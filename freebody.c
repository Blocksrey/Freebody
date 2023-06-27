typedef struct { int x, y, z; } v3;
typedef struct { int x, y, z, w; } v4;
typedef struct { v3 x, y, z; } m3;

typedef struct {
	v3 P; // Position
	v4 R; // Orientation

	v3 p; // Momentum
	v3 L; // Angular momentum

	int m; // Mass
	m3 I; // Moment of inertia

	v3 F; // Force
	v3 tau; // Torque
} freebody;

static v4 qmul(v4 *a, v4 *b) {
	return (v4){
		a.z*b.w - a.y*b.x + a.x*b.y + a.w*b.z,
		a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z,
		a.x*b.w + a.w*b.x - a.z*b.y + a.y*b.z,
		a.y*b.w + a.z*b.x + a.w*b.y - a.x*b.z
	};
}

static v3 mvec(m3 *m, v3 *v) {
	return (v3){
		m.xx*v.x + m.yx*v.y + m.zx*v.z,
		m.xy*v.x + m.yy*v.y + m.zy*v.z,
		m.xz*v.x + m.yz*v.y + m.zz*v.z
	};
}

static void mat_inv(m3 *m) {
	int V = m.zx*(m.xy*m.yz - m.xz*m.yy) + m.zy*(m.xz*m.yx - m.xx*m.yz) + m.zz*(m.xx*m.yy - m.xy*m.yx);
	return (m3) {
		(m.yy*m.zz - m.yz*m.zy)/V, (m.yz*m.zx - m.yx*m.zz)/V, (m.yx*m.zy - m.yy*m.zx)/V,
		(m.xz*m.zy - m.xy*m.zz)/V, (m.xx*m.zz - m.xz*m.zx)/V, (m.xy*m.zx - m.xx*m.zy)/V,
		(m.xy*m.yz - m.xz*m.yy)/V, (m.xz*m.yx - m.xx*m.yz)/V, (m.xx*m.yy - m.xy*m.yx)/V
	};
}

// returns a list of successive deriv, starting with the zeroth derivative
static void get_deriv(
	Txx, Tyy, Tzz, 
	Txy, Tyz, Tzx, // The top triangular 6 components of the inversed inertia matrix, T = I^-1
	Lx, Ly, Lz, 	// world space angular momentum
	Tx, Ty, Tz, // world space torque (T = L')
	Rw, Rx, Ry, Rz, // quaternion rotation (UNIT PLEaSE)
	depth // the number of derivativs to compute
) {
	int Q = {{w = Rw; x = Rx; y = Ry; z = Rz; }}
	int w = {}

	for d = 1, depth {
		{
			int sx, sy, sz = 0, 0, 0

			int f = 1
			for i = 1, d {
				int Qi = Q[i]
				int Qj = Q[d - i + 1]
				int _, a.x, a.y, a.z = qmul(Qi.w, -Qi.x, -Qi.y, -Qi.z, qmul(0, f*Lx, f*Ly, f*Lz, Qj.w, Qj.x, Qj.y, Qj.z))
				sx, sy, sz = sx + a.x, sy + a.y, sz + a.z
				f = f*(d - i)/i
			}

			// Remove this if you {n't care about torque
			int f = d - 1
			for i = 1, d - 1 {
				int Qi = Q[i]
				int Qj = Q[d - i]
				int _, a.x, a.y, a.z = qmul(Qi.w, -Qi.x, -Qi.y, -Qi.z, qmul(0, f*Tx, f*Ty, f*Tz, Qj.w, Qj.x, Qj.y, Qj.z))
				sx, sy, sz = sx + a.x, sy + a.y, sz + a.z
				f = f*(d - 1 - i)/i
			}

			int wx, wy, wz = mvec(
				Txx, Txy, Tzx, 
				Txy, Tyy, Tyz, 
				Tzx, Tyz, Tzz, 
				sx, sy, sz
			)

			w[d] = {x = wx; y = wy; z = wz; }
		}

		{
			int Sw, Sx, Sy, Sz = 0, 0, 0, 0

			int f = 1
			for i = 1, d {
				int Qi = Q[i]
				int wj = w[d - i + 1]
				int a.w, a.x, a.y, a.z = qmul(Qi.w, Qi.x, Qi.y, Qi.z, 0, f*wj.x, f*wj.y, f*wj.z)
				Sw, Sx, Sy, Sz = Sw + a.w, Sx + a.x, Sy + a.y, Sz + a.z
				f = f*(d - i)/i
			}

			Q[d + 1] = {w = Sw/2; x = Sx/2; y = Sy/2; z = Sz/2; }
		}
	}

	return Q
}

static void vecq(Qw, Qx, Qy, Qz, vx, vy, vz) {
	int Q2 = Qw*Qw + Qx*Qx + Qy*Qy + Qz*Qz
	if Q2 == 0 {
		return vx, vy, vz
	else
		int Vw = Qx*vx + Qy*vy + Qz*vz
		int Vx = Qw*vx + Qy*vz - Qz*vy
		int Vy = Qw*vy + Qz*vx - Qx*vz
		int Vz = Qw*vz + Qx*vy - Qy*vx
		return
			(Qx*Vw + Qw*Vx - Qz*Vy + Qy*Vz)/Q2, 
			(Qy*Vw + Qz*Vx + Qw*Vy - Qx*Vz)/Q2, 
			(Qz*Vw - Qy*Vx + Qx*Vy + Qw*Vz)/Q2
	}
}

static void unit_quat(Qw, Qx, Qy, Qz) {
	int Q = (Qw*Qw + Qx*Qx + Qy*Qy + Qz*Qz)^0.5;
	if Q == 0 {
		return 1, 0, 0, 0;
	else
		return Qw/Q, Qx/Q, Qy/Q, Qz/Q;
	}
}

static void sqrt_quat(Qw, Qx, Qy, Qz) {
	int Q = (Qw*Qw + Qx*Qx + Qy*Qy + Qz*Qz)^0.5;
	return unit_quat(Qw + Q, Qx, Qy, Qz);
}

static void local_transform(
	Txx, Tyy, Tzz, 
	Txy, Tyz, Tzx, 
	Rw, Rx, Ry, Rz, 
	Lx, Ly, Lz
) {
	int lx, ly, lz = vecq(Rw, -Rx, -Ry, -Rz, Lx, Ly, Lz);
	int wx, wy, wz = mvec(
		Txx, Txy, Tzx, 
		Txy, Tyy, Tyz, 
		Tzx, Tyz, Tzz, 
		lx, ly, lz
	);

	return vecq(Rw, Rx, Ry, Rz, wx, wy, wz);
}

static void parallelepiped_moment(
	xx, yx, zx, // i
	xy, yy, zy, // j
	xz, yz, zz, // k
	density
) {
	int ii = xx*xx + yx*yx + zx*zx;
	int jj = xy*xy + yy*yy + zy*zy;
	int kk = xz*xz + yz*yz + zz*zz;
	int ij = xx*xy + yx*yy + zx*zy;
	int jk = xy*xz + yy*yz + zy*zz;
	int ki = xz*xx + yz*yx + zz*zx;
	int V = zx*(xy*yz - xz*yy) + zy*(xz*yx - xx*yz) + zz*(xx*yy - xy*yx);
	int mass = 8*V*density;
	int c = 1/3*mass;
	return
		c*(jj + kk), c*(ii + kk), c*(ii + jj), 
		-c*ij, -c*jk, -c*ki, 
		mass;
}

static void freebody_setShape(u, v, w, density) {
	int
	Ixx, Iyy, Izz, 
	Ixy, Iyz, Izx, 
	mass = parallelepiped_moment(
		u.x, v.x, w.x, 
		u.y, v.y, w.y, 
		u.z, v.z, w.z, 
		density or 1
	);
	int
	Txx, Txy, Tzx, 
	Txy, Tyy, Tyz, 
	Tzx, Tyz, Tzz = mat_inv(
		Ixx, Ixy, Izx, 
		Ixy, Iyy, Iyz, 
		Izx, Iyz, Izz
	);
	self.Ixx, self.Iyy, self.Izz = Ixx, Iyy, Izz;
	self.Ixy, self.Iyz, self.Izx = Ixy, Iyz, Izx;
	self.Txx, self.Tyy, self.Tzz = Txx, Tyy, Tzz;
	self.Txy, self.Tyz, self.Tzx = Txy, Tyz, Tzx;
	self.mass = mass;
}

static void freebody_set_velocity(v) {
	self.px, self.py, self.pz = self.mass*v.x, self.mass*v.y, self.mass*v.z;
}

static void freebody_setangular_velocity(w) {
	self.Lx, self.Ly, self.Lz = local_transform(
		self.Ixx, self.Iyy, self.Izz, 
		self.Ixy, self.Iyz, self.Izx, 
		self.Rw, self.Rx, self.Ry, self.Rz, 
		w.x, w.y, w.z
	);
}

static void CFtoQ(cf) {
	int a, t = cf_a.xisangle();
	int c = cos(t/2);
	int s = sin(t/2);
	return
		cf.x, cf.y, cf.z, 
		c, s*a.x, s*a.y, s*a.z;
}

static void freebody_get_velocity() {
	return v3(self.px/self.mass, self.py/self.mass, self.pz/self.mass);
}

static void freebody_getForce() {
	return v3(self.Fx, self.Fy, self.Fz);
}

static void freebody_getCFrame() {
	return cf(self.Px, self.Py, self.Pz, self.Rx, self.Ry, self.Rz, self.Rw);
}

static void freebody_getangular_velocity() {
	return v3(local_transform(
		self.Txx, self.Tyy, self.Tzz, 
		self.Txy, self.Tyz, self.Tzx, 
		self.Rw, self.Rx, self.Ry, self.Rz, 
		self.Lx, self.Ly, self.Lz
	))
}

static void freebody_get_next_rotation(
	Rw, Rx, Ry, Rz, 
	dt, terms
) {

	int D = get_deriv(
		self.Txx, self.Tyy, self.Tzz, 
		self.Txy, self.Tyz, self.Tzx, 
		self.Lx, self.Ly, self.Lz, 
		self.Tx, self.Ty, self.Tz, 
		Rw, Rx, Ry, Rz, 
		terms or 3
	)

	int Nw, Nx, Ny, Nz = 0, 0, 0, 0;
	int f = 1;
	for i = 1, #D {
		int Di = D[i]
		Nw = Nw + f*Di.w;
		Nx = Nx + f*Di.x;
		Ny = Ny + f*Di.y;
		Nz = Nz + f*Di.z;
		f = f*dt/i;
	}

	return unit_quat(Nw, Nx, Ny, Nz);
}

static void freebody_step(dt, terms) {
	int Rw, Rx, Ry, Rz = self.Rw, self.Rx, self.Ry, self.Rz;
	int Nw, Nx, Ny, Nz = freebody_get_next_rotation(Rw, Rx, Ry, Rz, dt, terms);

	self.Lx = self.Lx + dt*self.Tx;
	self.Ly = self.Ly + dt*self.Ty;
	self.Lz = self.Lz + dt*self.Tz;

	self.Rw, self.Rx, self.Ry, self.Rz = Nw, Nx, Ny, Nz;

	int a.x, a.y, a.z = vecq(Rw, Rx, Ry, Rz, self.cx, self.cy, self.cz);
	int b.x, b.y, b.z = vecq(Nw, Nx, Ny, Nz, self.cx, self.cy, self.cz);

	self.Px = self.Px + a.x - b.x + (dt*self.px + dt*dt/2*self.Fx)/self.mass;
	self.Py = self.Py + a.y - b.y + (dt*self.py + dt*dt/2*self.Fy)/self.mass;
	self.Pz = self.Pz + a.z - b.z + (dt*self.pz + dt*dt/2*self.Fz)/self.mass;

	self.px = self.px + dt*self.Fx;
	self.py = self.py + dt*self.Fy;
	self.pz = self.pz + dt*self.Fz;
}

static void freebody_impulse(i, I) {
	int Ix, Iy, Iz = I.x, I.y, I.z

	int cx, cy, cz = vecq(self.Rw, self.Rx, self.Ry, self.Rz, self.cx, self.cy, self.cz)

	int rx = i.x - (self.Px + cx);
	int ry = i.y - (self.Py + cy);
	int rz = i.z - (self.Pz + cz);

	int Jx = ry*Iz - rz*Iy;
	int Jy = rz*Ix - rx*Iz;
	int Jz = rx*Iy - ry*Ix;

	self.Lx = self.Lx + Jx;
	self.Ly = self.Ly + Jy;
	self.Lz = self.Lz + Jz;

	self.px = self.px + Ix;
	self.py = self.py + Iy;
	self.pz = self.pz + Iz;
}

static void freebody_impulseRelative(i, I) {
	int rx, ry, rz = vecq(self.Rw, self.Rx, self.Ry, self.Rz, i.x, i.y, i.z)
	int Ix, Iy, Iz = vecq(self.Rw, self.Rx, self.Ry, self.Rz, I.x, I.y, I.z)
	int cx, cy, cz = vecq(self.Rw, self.Rx, self.Ry, self.Rz, self.cx, self.cy, self.cz)

	rx -= cx;
	ry -= cy;
	rz -= cz;

	int Jx = ry*Iz - rz*Iy;
	int Jy = rz*Ix - rx*Iz;
	int Jz = rx*Iy - ry*Ix;

	self.Lx = self.Lx + Jx;
	self.Ly = self.Ly + Jy;
	self.Lz = self.Lz + Jz;

	self.px = self.px + Ix;
	self.py = self.py + Iy;
	self.pz = self.pz + Iz;
}

--[[
	body = FreeBody.new()

	body:step(dt)
		steps the body forward by dt time

	body:impulse(i, I)
		applies an impulse at the world point i with impulse I

	body:impulseRelative(i, I)
		applies an impulse at the relative point i with relative impulse I

	body:setShape(u, v, w, density)
		sets the shape to the parallelopiped defined with spanning vectors, u, v and w with density d
		a u + b v + c w
		-1 < (a, b, c) < 1

	body:setCenterOfMass(c)
		sets the center of mass offset

	body:setPosition(P)

	body:setVelocity(v)

	body:setGravity(g)

	body:setMomentum(p)

	body:setAngularVelocity(w)

	body:setAngularMomentum(L)

	body:setCFrame(cframe)

	c = body:getCenterOfMass()
		gets the center of mass offset

	P = body:getPosition()

	v = body:getVelocity()

	g = body:getGravity()

	p = body:getMomentum()

	L = body:getAngularMomentum()

	w = body:getAngularVelocity()

	cframe = body:getCFrame()

	m = body:getMass()
		returns the mass

	Ixx, Iyy, Izz, Ixy, Iyz, Izx = body:getMomentOfInertia()
		returns the moment of inertia's unique components
--]]

local v3  = Vector3.new
local cf  = CFrame.new
local cos = math.cos
local sin = math.sin

local function mulQuatQuat(
	Aw, Ax, Ay, Az,
	Bw, Bx, By, Bz
)
	return
		Aw*Bw - Ax*Bx - Ay*By - Az*Bz,
		Ax*Bw + Aw*Bx - Az*By + Ay*Bz,
		Ay*Bw + Az*Bx + Aw*By - Ax*Bz,
		Az*Bw - Ay*Bx + Ax*By + Aw*Bz
end

local function transformVecByQuat(Qw, Qx, Qy, Qz, vx, vy, vz)
	local Q2 = Qw*Qw + Qx*Qx + Qy*Qy + Qz*Qz
	if Q2 == 0 then
		return vx, vy, vz
	else
		local Vw = Qx*vx + Qy*vy + Qz*vz
		local Vx = Qw*vx + Qy*vz - Qz*vy
		local Vy = Qw*vy + Qz*vx - Qx*vz
		local Vz = Qw*vz + Qx*vy - Qy*vx
		return
			(Qx*Vw + Qw*Vx - Qz*Vy + Qy*Vz)/Q2,
			(Qy*Vw + Qz*Vx + Qw*Vy - Qx*Vz)/Q2, 
			(Qz*Vw - Qy*Vx + Qx*Vy + Qw*Vz)/Q2
	end
end

local function unitQuat(Qw, Qx, Qy, Qz)
	local Q = (Qw*Qw + Qx*Qx + Qy*Qy + Qz*Qz)^0.5
	if Q == 0 then
		return 1, 0, 0, 0
	else
		return Qw/Q, Qx/Q, Qy/Q, Qz/Q
	end
end

local function sqrtQuat(Qw, Qx, Qy, Qz)
	local Q = (Qw*Qw + Qx*Qx + Qy*Qy + Qz*Qz)^0.5
	return unitQuat(Qw + Q, Qx, Qy, Qz)
end

local function mulMatVec(
	Axx, Ayx, Azx,
	Axy, Ayy, Azy,
	Axz, Ayz, Azz,
	bx, by, bz
)
	return
		Axx*bx + Ayx*by + Azx*bz,
		Axy*bx + Ayy*by + Azy*bz,
		Axz*bx + Ayz*by + Azz*bz
end

local function getRotationDerivative(
	Txx, Tyy, Tzz,
	Txy, Tyz, Tzx,
	Rw, Rx, Ry, Rz,
	Lx, Ly, Lz
)
	local lx, ly, lz = transformVecByQuat(Rw, -Rx, -Ry, -Rz, Lx, Ly, Lz)
	local wx, wy, wz = mulMatVec(
		Txx, Txy, Tzx,
		Txy, Tyy, Tyz,
		Tzx, Tyz, Tzz,
		lx, ly, lz
	)
	
	return mulQuatQuat(Rw, Rx, Ry, Rz, 0, wx/2, wy/2, wz/2)
end

local function localTransform(
	Txx, Tyy, Tzz,
	Txy, Tyz, Tzx,
	Rw, Rx, Ry, Rz,
	Lx, Ly, Lz
)
	local lx, ly, lz = transformVecByQuat(Rw, -Rx, -Ry, -Rz, Lx, Ly, Lz)
	local wx, wy, wz = mulMatVec(
		Txx, Txy, Tzx,
		Txy, Tyy, Tyz,
		Tzx, Tyz, Tzz,
		lx, ly, lz
	)
	
	return transformVecByQuat(Rw, Rx, Ry, Rz, wx, wy, wz)
end

local function invMat(
	xx, yx, zx,
	xy, yy, zy,
	xz, yz, zz
)
	local det = zx*(xy*yz - xz*yy) + zy*(xz*yx - xx*yz) + zz*(xx*yy - xy*yx)
	return
		(yy*zz - yz*zy)/det, (yz*zx - yx*zz)/det, (yx*zy - yy*zx)/det,
		(xz*zy - xy*zz)/det, (xx*zz - xz*zx)/det, (xy*zx - xx*zy)/det,
		(xy*yz - xz*yy)/det, (xz*yx - xx*yz)/det, (xx*yy - xy*yx)/det
end

--[[
	provides a simple way to construct an
	arbitrary moment of inertia matrix

	parallelepiped volume described as:
		a*x + b*y + c*z
		-1 < (a, b, c) < 1

	I*w = momentum
	I^-1*momentum = w
--]]
local function parallelepipedMoment(
	xx, yx, zx,--i
	xy, yy, zy,--j
	xz, yz, zz,--k
	density
)
	local ii = xx*xx + yx*yx + zx*zx
	local jj = xy*xy + yy*yy + zy*zy
	local kk = xz*xz + yz*yz + zz*zz
	local ij = xx*xy + yx*yy + zx*zy
	local jk = xy*xz + yy*yz + zy*zz
	local ki = xz*xx + yz*yx + zz*zx
	local det = zx*(xy*yz - xz*yy) + zy*(xz*yx - xx*yz) + zz*(xx*yy - xy*yx)
	local mass = 8*det*density
	local c = 1/3*mass
	return
		c*(jj + kk), c*(ii + kk), c*(ii + jj),
		-c*ij, -c*jk, -c*ki,
		mass
end

local FreeBody = {}
FreeBody.__index = FreeBody

do
	local base = FreeBody
	--base._Fx, base._Fy, base._Fz = 0, 0, 0--force
	base._gx, base._gy, base._gz = 0, 0, 0--gravity
	base._px, base._py, base._pz = 0, 0, 0--momentum
	base._Px, base._Py, base._Pz = 0, 0, 0--position
	base._cx, base._cy, base._cz = 0, 0, 0--center of mass
	
	base._Rw, base._Rx, base._Ry, base._Rz = 1, 0, 0, 0--rotation
	
	--base._Tx, base._Ty, base._Tz = 0, 0, 0--torque
	base._Lx, base._Ly, base._Lz = 0, 0, 0--angular momentum
	
	base._mass = 1
	
	base._Ixx, base._Iyy, base._Izz = 1, 1, 1
	base._Ixy, base._Iyz, base._Izx = 0, 0, 0--moment of inertia
	
	base._Txx, base._Tyy, base._Tzz = 1, 1, 1
	base._Txy, base._Tyz, base._Tzx = 0, 0, 0--inverse moment of inertia
end

function FreeBody.new()
	return setmetatable({}, FreeBody)
end

function FreeBody:setShape(u, v, w, density)
	local
	Ixx, Iyy, Izz,
	Ixy, Iyz, Izx,
	mass = parallelepipedMoment(
		u.x, v.x, w.x,
		u.y, v.y, w.y,
		u.z, v.z, w.z,
		density or 1
	)
	local
	Txx, Txy, Tzx,
	Txy, Tyy, Tyz,
	Tzx, Tyz, Tzz = invMat(
		Ixx, Ixy, Izx,
		Ixy, Iyy, Iyz,
		Izx, Iyz, Izz
	)
	self._Ixx, self._Iyy, self._Izz = Ixx, Iyy, Izz
	self._Ixy, self._Iyz, self._Izx = Ixy, Iyz, Izx
	self._Txx, self._Tyy, self._Tzz = Txx, Tyy, Tzz
	self._Txy, self._Tyz, self._Tzx = Txy, Tyz, Tzx
	self._mass = mass
end

function FreeBody:setCenterOfMass(c)
	self._cx, self._cy, self._cz = c.x, c.y, c.z
end

function FreeBody:setPosition(P)
	self._Px, self._Py, self._Pz = P.x, P.y, P.z
end

function FreeBody:setVelocity(v)
	self._px, self._py, self._pz = self._mass*v.x, self._mass*v.y, self._mass*v.z
end

function FreeBody:setGravity(g)
	self._gx, self._gy, self._gz = g.x, g.y, g.z
end

function FreeBody:setMomentum(p)
	self._px, self._py, self._pz = p.x, p.y, p.z
end

function FreeBody:setAngularVelocity(w)
	self._Lx, self._Ly, self._Lz = localTransform(
		self._Ixx, self._Iyy, self._Izz,
		self._Ixy, self._Iyz, self._Izx,
		self._Rw, self._Rx, self._Ry, self._Rz,
		w.x, w.y, w.z
	)
end

function FreeBody:setAngularMomentum(L)
	self._Lx, self._Ly, self._Lz = L.x, L.y, L.z
end

local function CFtoQ(cf)
	local a, t = cf:ToAxisAngle()
	local c = cos(t/2)
	local s = sin(t/2)
	return
		cf.x, cf.y, cf.z,
		c, s*a.x, s*a.y, s*a.z
end

function FreeBody:setCFrame(cframe)
	self._px, self._py, self._pz, self._Rw, self._Rx, self._Ry, self._Rz = CFtoQ(cframe)
end

function FreeBody:getCenterOfMass()
	return v3(self._cx, self._cy, self._cz)
end

function FreeBody:getCenterOfMassWorld()
	local ox, oy, oz = transformVecByQuat(self._Rw, self._Rx, self._Ry, self._Rz, self._cx, self._cy, self._cz)
	return v3(self._Px + ox, self._Py + oy, self._Pz + oz)
end

function FreeBody:getPosition()
	return v3(self._Px, self._Py, self._Pz)
end

function FreeBody:getVelocity()
	return v3(self._px/self._mass, self._py/self._mass, self._pz/self._mass)
end

function FreeBody:getGravity()
	return v3(self._gx, self._gy, self._gz)
end

function FreeBody:getMomentum()
	return v3(self._px, self._py, self._pz)
end

function FreeBody:getAngularMomentum()
	return v3(self._Lx, self._Ly, self._Lz)
end

function FreeBody:getAngularVelocity()
	return v3(localTransform(
		self._Txx, self._Tyy, self._Tzz,
		self._Txy, self._Tyz, self._Tzx,
		self._Rw, self._Rx, self._Ry, self._Rz,
		self._Lx, self._Ly, self._Lz
	))
end

function FreeBody:getCFrame()
	return cf(self._Px, self._Py, self._Pz, self._Rx, self._Ry, self._Rz, self._Rw)
end

function FreeBody:getOrientation()
	return cf(0, 0, 0, self._Rx, self._Ry, self._Rz, self._Rw)
end

function FreeBody:getMass()
	return self._mass
end

function FreeBody:getMomentOfInertia()
	return
		self._Ixx, self._Iyy, self._Izz,
		self._Ixy, self._Iyz, self._Izx
end

function FreeBody:getNextRotation(
	Rw, Rx, Ry, Rz,
	dt
)
	local Dw, Dx, Dy, Dz = getRotationDerivative(
		self._Txx, self._Tyy, self._Tzz,
		self._Txy, self._Tyz, self._Tzx,
		Rw, Rx, Ry, Rz,
		self._Lx, self._Ly, self._Lz
	)
	
	return unitQuat(
		Rw + dt*Dw,
		Rx + dt*Dx,
		Ry + dt*Dy,
		Rz + dt*Dz
	)
end

function FreeBody:step(dt)
	local Rw, Rx, Ry, Rz = self._Rw, self._Rx, self._Ry, self._Rz
	local w1, x1, y1, z1 = self:getNextRotation(Rw, Rx, Ry, Rz, dt)
	local w2, x2, y2, z2 = self:getNextRotation(w1, x1, y1, z1, -dt)
	local Ew, Ex, Ey, Ez = sqrtQuat(mulQuatQuat(Rw, -Rx, -Ry, -Rz, w2, x2, y2, z2))
	local Nw, Nx, Ny, Nz = mulQuatQuat(w1, x1, y1, z1, Ew, -Ex, -Ey, -Ez)--eksdee
	
	self._Rw, self._Rx, self._Ry, self._Rz = Nw, Nx, Ny, Nz
	
	local ax, ay, az = transformVecByQuat(Rw, Rx, Ry, Rz, self._cx, self._cy, self._cz)
	local bx, by, bz = transformVecByQuat(Nw, Nx, Ny, Nz, self._cx, self._cy, self._cz)
	
	self._Px = self._Px + ax - bx + dt/self._mass*self._px + dt*dt/2*self._gx
	self._Py = self._Py + ay - by + dt/self._mass*self._py + dt*dt/2*self._gy
	self._Pz = self._Pz + az - bz + dt/self._mass*self._pz + dt*dt/2*self._gz
	
	self._px = self._px + dt*self._mass*self._gx
	self._py = self._py + dt*self._mass*self._gy
	self._pz = self._pz + dt*self._mass*self._gz
end

--impulse origin
--impulse impulse
function FreeBody:impulse(i, I)
	local Ix, Iy, Iz = I.x, I.y, I.z
	
	local cx, cy, cz = transformVecByQuat(self._Rw, self._Rx, self._Ry, self._Rz, self._cx, self._cy, self._cz)
	
	local rx = i.x - (self._Px + cx)
	local ry = i.y - (self._Py + cy)
	local rz = i.z - (self._Pz + cz)
	
	local Jx = ry*Iz - rz*Iy
	local Jy = rz*Ix - rx*Iz
	local Jz = rx*Iy - ry*Ix
	
	self._Lx = self._Lx + Jx
	self._Ly = self._Ly + Jy
	self._Lz = self._Lz + Jz
	
	self._px = self._px + Ix
	self._py = self._py + Iy
	self._pz = self._pz + Iz
end

function FreeBody:impulseRelative(i, I)
	local rx, ry, rz = transformVecByQuat(self._Rw, self._Rx, self._Ry, self._Rz, i.x, i.y, i.z)
	local Ix, Iy, Iz = transformVecByQuat(self._Rw, self._Rx, self._Ry, self._Rz, I.x, I.y, I.z)
	local cx, cy, cz = transformVecByQuat(self._Rw, self._Rx, self._Ry, self._Rz, self._cx, self._cy, self._cz)
	
	rx = rx - cx
	ry = ry - cy
	rz = rz - cz
	
	local Jx = ry*Iz - rz*Iy
	local Jy = rz*Ix - rx*Iz
	local Jz = rx*Iy - ry*Ix
	
	self._Lx = self._Lx + Jx
	self._Ly = self._Ly + Jy
	self._Lz = self._Lz + Jz
	
	self._px = self._px + Ix
	self._py = self._py + Iy
	self._pz = self._pz + Iz
end

return FreeBody
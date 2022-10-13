--Created by Trey Reynolds (AxisAngle) & modified by Jeff Skinner (Blocksrey)

--Transitional/Rotational
--Position: P/R
--Momentum: p/L
--Mass: m/I
--Velocity: v/omega
--Force: F/tau
--Impulse: J/H
--Impulse: delta_p/delta_L

local vec3 = require('vec3')
local mat3 = require('mat3')
local quat = require('quat')

local freebody = {}
freebody.__index = freebody

function freebody.new()
	local self = setmetatable({}, freebody)

	self.P = vec3.new(0, 0, 0)--Position
	self.R = quat.new(1, 0, 0, 0)--Orientation

	self.p = vec3.new(0, 0, 0)--Momentum
	self.L = vec3.new(0, 0, 0)--Angular momentum

	self.m = 1--Mass
	self:set_I(mat3.identity)

	self.F = vec3.new(0, 0, 0)--Force
	self.tau = vec3.new(0, 0, 0)--Torque

	return self
end

--Set shape to an arbitrary parallelepiped
function freebody:set_shape(M, D)
	local m = 8*D*M:det()--8*density*volume
	local c = 1/3*m

	local i = M*vec3.new(1, 0, 0)
	local j = M*vec3.new(0, 1, 0)
	local k = M*vec3.new(0, 0, 1)

	local ii = i:dot(i)
	local jj = j:dot(j)
	local kk = k:dot(k)

	local ij = i:dot(j)
	local jk = j:dot(k)
	local ki = k:dot(i)

	local I = mat3.new(
		c*(jj + kk), -c*ij, -c*ki,
		-c*ij, c*(ii + kk), -c*jk,
		-c*ki, -c*jk, c*(ii + jj)
	)

	self.m = m
	self:set_I(I)
end

function freebody:set_omega(omega)
	self.L = self.R*(self.I*(self.R:inv()*omega))
end

function freebody:set_I(I)
	self.I = I
	self.i = I:inv()
end

local function rotation_derivatives(i, L, tau, R, n_terms)
	local Q = {R}
	local w = {}

	for d = 1, n_terms do
		local s = vec3.new(0, 0, 0)

		---[[
		local f = 1
		for i = 1, d do
			s = s + Q[i]:inv()*quat.new(0, (f*L):dump())*Q[d - i + 1]
			f = f*(d - i)/i
		end
		--]]

		---[[
		local f = d - 1
		for i = 1, d - 1 do
			s = s + Q[i]:inv()*quat.new(0, (f*tau):dump())*Q[d - i]
			f = f*(d - i - 1)/i
		end
		--]]

		w[d] = i*s

		local S = quat.new(0, 0, 0, 0)

		---[[
		local f = 1
		for i = 1, d do
			S = S + Q[i]*quat.new(0, (f*w[d - i + 1]):dump())
			f = f*(d - i)/i
		end
		--]]

		Q[d + 1] = 0.5*S
	end

	return Q
end

function freebody:step(dt, n_terms)
	local D = rotation_derivatives(self.i, self.L, self.tau, self.R, n_terms)
	local N = quat.new(0, 0, 0, 0)
	local f = 1
	for i = 1, #D do
		N = N + f*D[i]
		f = f*dt/i
	end
	self.R = N:unit()
	self.L = self.L + dt*self.tau
	self.P = self.P + dt/self.m*self.p
end

function freebody:impulse(P, J)
	self.p = self.p + J
	self.L = self.L + (P - self.P):cross(J)
end

function freebody:angular_impulse(H)
	self.L = self.L + H
end

return freebody
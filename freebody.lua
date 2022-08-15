--Transitional/Rotational
--Position: P/R
--Momentum: p/L
--Mass: m/I
--Velocity: v/omega
--Force: F/tau
--Impulse: J/H
--Impulse: delta_p/delta_L

local Vec3 = require(script.Parent.Vec3)
local Mat3 = require(script.Parent.Mat3)
local Quat = require(script.Parent.Quat)

local FreeBody = {}
FreeBody.__index = FreeBody

function FreeBody.new()
	local self = setmetatable({}, FreeBody)
	
	self._P = Vec3.new()--Position
	self._R = Quat.new()--Orientation
	
	self._p = Vec3.new()--Momentum
	self._L = Vec3.new()--Angular momentum
	
	self._m = 1--Mass
	self:SetI(Mat3.identity)
	
	self._F = Vec3.new()--Force
	self._tau = Vec3.new()--Torque
	
	return self
end

function FreeBody:GetP()
	return self._P
end

function FreeBody:SetP(P)
	self._P = P
end

function FreeBody:GetR()
	return self._R
end

function FreeBody:SetR(R)
	self._R = R
end

function FreeBody:Getp()
	return self._p
end

function FreeBody:Setp(p)
	self._p = p
end

function FreeBody:GetL()
	return self._L
end

function FreeBody:SetL(L)
	self._L = L
end

function FreeBody:SetShape(M, D)
	local det = M:Det()
	local m = 8*det*D
	local c = 1/3*m
	
	local i = M*Vec3.new(1, 0, 0)
	local j = M*Vec3.new(0, 1, 0)
	local k = M*Vec3.new(0, 0, 1)

	local ii = i:Dot(i)
	local jj = j:Dot(j)
	local kk = k:Dot(k)

	local ij = i:Dot(j)
	local jk = j:Dot(k)
	local ki = k:Dot(i)
	
	local I = Mat3.new(
		c*(jj + kk), -c*ij, -c*ki,
		-c*ij, c*(ii + kk), -c*jk,
		-c*ki, -c*jk, c*(ii + jj)
	)
	
	self._m = m
	self:SetI(I)
end

local function localTransform(T, R, L)
	return R*(T*(R:Inverse()*L))
end

function FreeBody:SetOmega(omega)
	self._L = localTransform(
		self._I,
		self._R,
		omega
	)
end

function FreeBody:SetI(I)
	self._I = I
	self._i = I:Inverse()
end

local function computeRotationDerivatives(i, L, tau, R, nTerms)
	local Q = {R}
	local w = {}
	
	for d = 1, nTerms do
		local s = Vec3.new()

		---[[
		local f = 1
		for i = 1, d do
			s = s + Q[i]:Inverse()*Quat.new(0, (f*L):Dump())*Q[d - i + 1]
			f = f*(d - i)/i
		end
		--]]

		---[[
		local f = d - 1
		for i = 1, d - 1 do
			s = s + Q[i]:Inverse()*Quat.new(0, (f*tau):Dump())*Q[d - i]
			f = f*(d - i - 1)/i
		end
		--]]
		
		w[d] = i*s
		
		local S = Quat.new(0)

		---[[
		local f = 1
		for i = 1, d do
			S = S + Q[i]*Quat.new(0, (f*w[d - i + 1]):Dump())
			f = f*(d - i)/i
		end
		--]]
		
		Q[d + 1] = 0.5*S
	end

	return Q
end

function FreeBody:Step(delta_t, nTerms)
	local D = computeRotationDerivatives(self._i, self._L, self._tau, self._R, nTerms)
	local N = Quat.new(0)
	local f = 1
	for i = 1, #D do
		N = N + f*D[i]
		f = f*delta_t/i
	end
	self._R = N:Unitize()
	self._L = self._L + delta_t*self._tau
	self._P = self._P + delta_t/self._m*self._p
end

function FreeBody:Impulse(P, J)
	self._p = self._p + J
	self._L = self._L + (P - self._P):Cross(J)
end

function FreeBody:AngularImpulse(H)
	self._L = self._L + H
end

return FreeBody

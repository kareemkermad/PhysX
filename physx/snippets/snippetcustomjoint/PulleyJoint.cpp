// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of NVIDIA CORPORATION nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
// OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Copyright (c) 2008-2024 NVIDIA Corporation. All rights reserved.
// Copyright (c) 2004-2008 AGEIA Technologies, Inc. All rights reserved.
// Copyright (c) 2001-2004 NovodeX AG. All rights reserved.  

#include "PulleyJoint.h"
#include "PxConstraint.h"
#include <assert.h>

using namespace physx;

//TAG:solverprepshader
static PxU32 solverPrep(Px1DConstraint* constraints,
						PxVec3p& body0WorldOffset,
						PxU32 maxConstraints,
						PxConstraintInvMassScale& invMassScale,
						const void* constantBlock,
						const PxTransform& bA2w,
						const PxTransform& bB2w,
						bool /*useExtendedLimits*/,
						PxVec3p& cA2wOut, PxVec3p& cB2wOut)
{		
	PX_UNUSED(maxConstraints);

	const PulleyJoint::PathJointData& data = *reinterpret_cast<const PulleyJoint::PathJointData*>(constantBlock);
	const curve* path = data.path;

	PxTransform32 cA2w, cB2w;
	physx::Ext::joint::ConstraintHelper ch(constraints, invMassScale, cA2w, cB2w, body0WorldOffset, data, bA2w, bB2w);

	physx::Ext::joint::applyNeighborhoodOperator(cA2w, cB2w);

	const float s = path->projected_length(bA2w.transformInv(cB2w.p));
	path->frame(s, cA2w);
	cA2w = bA2w.transform(cA2w);

	/*
	const bool closed = path->closed();
	bool limitEnabled = data.joint_flags & PxPathJointFlag::eLIMIT_ENABLED;
	const PxJointLinearLimitPair& limit = data.limit;
	float upper = limit.upper + data.zero_position; // Map limits to the interval [0, L]
	float lower = limit.lower + data.zero_position;

	if (closed == false)
	{
		if (limitEnabled)
		{
			upper = fminf(upper, path->length());
			lower = fmaxf(lower, 0.0F);
		}
		else
		{
			upper = path->length();
			lower = 0.0F;
			limitEnabled = true;
		}
	}

	const bool limitIsLocked = limitEnabled && lower >= upper;
	*/

	const bool limitIsLocked = false;

	const PxVec3 bOriginInA = cA2w.transformInv(cB2w.p);

	const PxU32 angle_constraints = PxU32(data.joint_flags & (PxPathJointFlag::eTANGENT_ANGLE_CONSTRAINT_ENABLED | PxPathJointFlag::eNORMAL_ANGLE_CONSTRAINT_ENABLED | PxPathJointFlag::eBITANGENT_ANGLE_CONSTRAINT_ENABLED)) >> 3;

	PxVec3 ra, rb, axis;
	ch.prepareLockedAxes(cA2w.q, cB2w.q, bOriginInA, limitIsLocked ? 7ul : 6ul, angle_constraints, ra, rb, &axis);
	cA2wOut = ra + bA2w.p;
	cB2wOut = rb + bB2w.p;

	if ((data.joint_flags & PxPathJointFlag::eDRIVE_ENABLED) == PxPathJointFlag::eDRIVE_ENABLED)
	{
		ch.linear(axis, -data.drive_velocity, data.drive_relative_position, data.drive_settings);
	}

	// Limits
	// The limits are linear!
	// We need to set upper and lower relative to s.
	// For example, if we are on a circle of circumference 1, we start at s=0 and our limits are lower=-0.1 and upper = 0.1.
	// Suppose we move in the negative direction and we are now at s=0.99. We can set lower=0.9 and upper = 1.1.

	/*
	if (limitEnabled && !limitIsLocked)
	{
		ch.linearLimit(axis, s, upper, limit);
		ch.linearLimit(-axis, -s, -lower, limit);
	}
	*/

	return ch.getCount();
}

static void visualize(PxConstraintVisualizer& viz,
					  const void* constantBlock,
					  const PxTransform& body0Transform,
					  const PxTransform& body1Transform,
					  PxU32	flags)
{
	PX_UNUSED(flags);

	const PulleyJoint::PathJointData& data = *reinterpret_cast<const PulleyJoint::PathJointData*>(constantBlock);
	PxTransform cA2w = body0Transform * data.c2b[0];
	PxTransform cB2w = body1Transform * data.c2b[1];
	viz.visualizeJointFrames(cA2w, cB2w);
	viz.visualizeJointFrames(PxTransform(PxIdentity), PxTransform(PxIdentity));
}

static PxConstraintShaderTable sShaderTable = { solverPrep, visualize, PxConstraintFlag::Enum(0) };

PxConstraintSolverPrep PulleyJoint::getPrep() const { return solverPrep; }

PulleyJoint::PulleyJoint(PxPhysics& physics, const physx::curve* path, PxRigidBody* body0, const PxTransform& localFrame0, PxRigidBody* body1, const PxTransform& localFrame1)
: PxJoint(PxCustomJointConcreteType::ePATH, PxBaseFlag::eOWNS_MEMORY | PxBaseFlag::eIS_RELEASABLE),
  m_data(path, PxJointLinearLimitPair(PxTolerancesScale()))
{
	m_constraint = physics.createConstraint(body1, body0, *this, sShaderTable, sizeof(PathJointData));

	m_bodies[0] = body1;
	m_bodies[1] = body0;

	// keep these around in case the CoM gets relocated
	m_local_poses[0] = localFrame1.getNormalized();
	m_local_poses[1] = localFrame0.getNormalized();

	// the data which will be fed to the joint solver and projection shaders
	m_data.invMassScale.linear0 = 1.0f;
	m_data.invMassScale.angular0 = 1.0f;
	m_data.invMassScale.linear1 = 1.0f;
	m_data.invMassScale.angular1 = 1.0f;
	m_data.c2b[0] = (body1 != nullptr) ? body1->getCMassLocalPose().transformInv(m_local_poses[0]) : m_local_poses[0];
	m_data.c2b[1] = (body0 != nullptr) ? body0->getCMassLocalPose().transformInv(m_local_poses[1]) : m_local_poses[1];
}

void PulleyJoint::release()
{
	m_constraint->release();
}

void PulleyJoint::setFlags(physx::PxPathJointFlags flags)
{
	m_data.joint_flags = flags;
	m_constraint->markDirty();
}

void PulleyJoint::setLimit(float lower, float upper)
{
	m_data.limit.lower = physx::PxReal(lower);
	m_data.limit.upper = physx::PxReal(upper);
	m_constraint->markDirty();
}

void PulleyJoint::setDriveRelativePosition(float relative_position)
{
	m_data.drive_relative_position = relative_position;
	m_constraint->markDirty();
}

void PulleyJoint::setDriveVelocity(float velocity)
{
	m_data.drive_velocity = velocity;
	m_constraint->markDirty();
}

void PulleyJoint::setDriveSettings(const physx::PxD6JointDrive& settings)
{
	m_data.drive_settings = settings;
	m_constraint->markDirty();
}

float PulleyJoint::getCurrentPosition()
{
	const physx::PxRigidBody* body0 = m_bodies[1];
	const physx::PxRigidBody* body1 = m_bodies[0];
	const physx::PxVec3 world_anchor = (body0 != nullptr) ? body0->getGlobalPose().transform(m_local_poses[1].p) : m_local_poses[1].p;
	const physx::PxVec3 local_anchor = (body1 != nullptr) ? body1->getGlobalPose().transformInv(world_anchor) : world_anchor;
	return m_data.path->projected_length(local_anchor);
}

const char* PulleyJoint::getConcreteTypeName() const
{
	return nullptr;
}

PxScene* PulleyJoint::getScene() const
{
	return nullptr;
}

physx::PxConstraint* PulleyJoint::getConstraint() const
{
	return m_constraint;
}

physx::PxConstraintFlags PulleyJoint::getConstraintFlags() const
{
	return m_constraint->getFlags();
}

void PulleyJoint::setConstraintFlags(PxConstraintFlags flags)
{
	m_constraint->setFlags(flags);
}

void PulleyJoint::setConstraintFlag(PxConstraintFlag::Enum flag, bool value)
{
	m_constraint->setFlag(flag, value);
}

const char* PulleyJoint::getName() const
{
	return nullptr;
}

void PulleyJoint::setName(const char* name)
{
	PX_UNUSED(name);
	return;
}

PxTransform PulleyJoint::getLocalPose(PxJointActorIndex::Enum actor) const
{
	switch (actor)
	{
	case PxJointActorIndex::eACTOR0:
		return m_local_poses[0];
	case PxJointActorIndex::eACTOR1:
		return m_local_poses[1];
	case PxJointActorIndex::COUNT:
		return PxTransform(PxIdentity);
	}

	return PxTransform(PxIdentity);
}

void PulleyJoint::setLocalPose(PxJointActorIndex::Enum actor, const PxTransform& localPose)
{
	switch (actor)
	{
	case PxJointActorIndex::eACTOR0:
		m_local_poses[0] = localPose;
		break;
	case PxJointActorIndex::eACTOR1:
		m_local_poses[1] = localPose;
		break;
	}
}

void PulleyJoint::setActors(PxRigidActor* actor0, PxRigidActor* actor1)
{
	PX_UNUSED(actor0);
	PX_UNUSED(actor1);
	return;
}

void PulleyJoint::getActors(PxRigidActor*& actor0, PxRigidActor*& actor1) const
{
	PX_UNUSED(actor0);
	PX_UNUSED(actor1);
	return;
}

PxTransform PulleyJoint::getRelativeTransform() const
{
	return PxTransform();
}

PxVec3 PulleyJoint::getRelativeLinearVelocity() const
{
	return PxVec3();
}

PxVec3 PulleyJoint::getRelativeAngularVelocity() const
{
	return PxVec3();
}

void PulleyJoint::getBreakForce(PxReal& force, PxReal& torque) const
{
	PX_UNUSED(force);
	PX_UNUSED(torque);
	return;
}

void PulleyJoint::setBreakForce(PxReal force, PxReal torque)
{
	PX_UNUSED(force);
	PX_UNUSED(torque);
	return;
}

void PulleyJoint::setInvMassScale0(PxReal invMassScale)
{
	PX_UNUSED(invMassScale);
	return;
}

PxReal PulleyJoint::getInvMassScale0() const
{
	return PxReal();
}

void PulleyJoint::setInvInertiaScale0(PxReal invInertiaScale)
{
	PX_UNUSED(invInertiaScale);
	return;
}

PxReal PulleyJoint::getInvInertiaScale0() const
{
	return PxReal();
}

void PulleyJoint::setInvMassScale1(PxReal invMassScale)
{
	PX_UNUSED(invMassScale);
	return;
}

PxReal PulleyJoint::getInvMassScale1() const
{
	return PxReal();
}

void PulleyJoint::setInvInertiaScale1(PxReal invInertiaScale)
{
	PX_UNUSED(invInertiaScale);
	return;
}

PxReal PulleyJoint::getInvInertiaScale1() const
{
	return PxReal();
}

void* PulleyJoint::prepareData()
{
	return &m_data;
}

void PulleyJoint::onConstraintRelease()
{
	delete this;
}

void PulleyJoint::onComShift(PxU32 actor)
{
	m_data.c2b[actor] = m_bodies[actor]->getCMassLocalPose().transformInv(m_local_poses[actor]); 
	m_constraint->markDirty();
}

void PulleyJoint::onOriginShift(const PxVec3& shift)
{
	PX_UNUSED(shift);
}

void* PulleyJoint::getExternalReference(PxU32& typeID)
{
	typeID = TYPE_ID;
	return this;
}


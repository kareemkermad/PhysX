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

	PxVec3 local_anchor_b = bA2w.transformInv(cB2w.p);
	const float s = path->projected_length(local_anchor_b);
	path->frame(s, cA2w);
	assert(cA2w.isSane());
	cA2w = bA2w.transform(cA2w);

	const PxVec3 bOriginInA = cA2w.transformInv(cB2w.p);

	const PxU32 angle_constraints = PxU32(data.joint_flags & (PxPathJointFlag::eTANGENT_ANGLE_CONSTRAINT_ENABLED | PxPathJointFlag::eNORMAL_ANGLE_CONSTRAINT_ENABLED | PxPathJointFlag::eBITANGENT_ANGLE_CONSTRAINT_ENABLED)) >> 3;

	PxVec3 ra, rb, axis;
	ch.prepareLockedAxes(cA2w.q, cB2w.q, bOriginInA, 6ul, angle_constraints, ra, rb, &axis);
	cA2wOut = ra + bA2w.p;
	cB2wOut = rb + bB2w.p;

	if ((data.joint_flags & PxPathJointFlag::eDRIVE_ENABLED) == PxPathJointFlag::eDRIVE_ENABLED)
	{
		const float drive_position_error = data.drive_position - s;
		ch.linear(axis, data.drive_velocity, drive_position_error, data.drive_settings);
	}

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

PulleyJoint::PulleyJoint(PxPhysics& physics, const physx::curve* path, PxRigidBody& body0, const PxTransform& localFrame0, PxRigidBody& body1, const PxTransform& localFrame1)
: m_data(path, PxJointLinearLimitPair(PxTolerancesScale()))
{
	m_constraint = physics.createConstraint(&body0, &body1, *this, sShaderTable, sizeof(PathJointData));

	m_bodies[0] = &body0;
	m_bodies[1] = &body1;

	// keep these around in case the CoM gets relocated
	m_local_poses[0] = localFrame0.getNormalized();
	m_local_poses[1] = localFrame1.getNormalized();

	// the data which will be fed to the joint solver and projection shaders
	m_data.invMassScale.linear0 = 1.0f;
	m_data.invMassScale.angular0 = 1.0f;
	m_data.invMassScale.linear1 = 1.0f;
	m_data.invMassScale.angular1 = 1.0f;
	m_data.c2b[0] = body0.getCMassLocalPose().transformInv(m_local_poses[0]);
	m_data.c2b[1] = body1.getCMassLocalPose().transformInv(m_local_poses[1]);
}

void PulleyJoint::release()
{
	m_constraint->release();
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


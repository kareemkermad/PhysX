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

#ifndef PULLEY_JOINT_H
#define PULLEY_JOINT_H

#include "PxPhysicsAPI.h"
#include "ExtConstraintHelper.h"

// a pulley joint constrains two actors such that the sum of their distances from their respective anchor points at their attachment points 
// is a fixed value (the parameter 'distance'). Only dynamic actors are supported.
//
// The constraint equation is as follows:
//
// |anchor0 - attachment0| + |anchor1 - attachment1| * ratio = distance
// 
// where 'ratio' provides mechanical advantage.
//
// The above equation results in a singularity when the anchor point is coincident with the attachment point; for simplicity
// the constraint does not attempt to handle this case robustly.

namespace physx
{
    struct PxPathJointFlag
    {
        enum Enum
        {
            eLIMIT_ENABLED = 1 << 1,
            eDRIVE_ENABLED = 1 << 2,
            eTANGENT_ANGLE_CONSTRAINT_ENABLED = 1 << 3,
            eNORMAL_ANGLE_CONSTRAINT_ENABLED = 1 << 4,
            eBITANGENT_ANGLE_CONSTRAINT_ENABLED = 1 << 5
        };
    };

    typedef PxFlags<PxPathJointFlag::Enum, PxU16> PxPathJointFlags;
    PX_FLAGS_OPERATORS(PxPathJointFlag::Enum, PxU16)

    class plane
    {
    public:
        inline plane() {}

        inline plane(const PxVec3& origin, const PxVec3& first_axis, const PxVec3& second_axis, const PxVec3& normal)
        : m_origin(origin),
          m_first_axis(first_axis),
          m_second_axis(second_axis),
          m_normal(normal)
        {
            // Nothing to do.
        }

    public:
        inline const PxVec3& origin() const
        {
            return m_origin;
        }

        inline const PxVec3& normal() const
        {
            return m_normal;
        }

        inline const PxVec3& first_axis() const
        {
            return m_first_axis;
        }

        inline const PxVec3& second_axis() const
        {
            return m_second_axis;
        }

    public:
        static inline float dot(PxVec3 a, PxVec3 b)
        {
            return a.x * b.x + a.y * b.y + a.z * b.z;
        }

        static inline float signed_distance(const plane& plane, const PxVec3& point)
        {
            return dot(point - plane.m_origin, plane.m_normal);
        }

        static inline PxVec3 project(const plane& plane, const PxVec3& point)
        {
            const PxVec3 v = point - plane.m_origin;
            const float px = dot(v, plane.m_first_axis);
            const float py = dot(v, plane.m_second_axis);
            return plane.m_origin + plane.m_first_axis * px + plane.m_second_axis * py;
        }

        static inline float projected_angle(const plane& plane, const PxVec3& point)
        {
            const PxVec3 o = plane.origin();
            const PxVec3 n = plane.normal();
            const PxVec3 f = plane.first_axis();
            const float tox = point.x - o.x;
            const float toy = point.y - o.y;
            const float toz = point.z - o.z;
            const float l = PxSqrt(tox * tox + toy * toy + toz * toz);
            const float costheta = (f.x * tox + f.y * toy + f.z * toz) / l;
            const float a = PxAcos(PxClamp(costheta, -1.0f, 1.0f));
            const float cx = f.y * toz - f.z * toy;
            const float cy = f.z * tox - f.x * toz;
            const float cz = f.x * toy - f.y * tox;
            const float d = n.x * cx + n.y * cy + n.z * cz;
            return (d >= 0) ? a : PxTwoPi - a;
        }

        static inline PxVec3 planar_coordinates(const plane& plane, const PxVec3& point)
        {
            const PxVec3 v = point - plane.m_origin;
            return PxVec3(dot(v, plane.m_first_axis), dot(v, plane.m_second_axis), dot(v, plane.m_normal));
        }

        static inline PxVec2 planar_coordinates_2d(const plane& plane, const PxVec3& point)
        {
            const PxVec3 v = point - plane.m_origin;
            return PxVec2(dot(v, plane.m_first_axis), dot(v, plane.m_second_axis));
        }

        static inline PxVec3 cartesian_coordinates(const plane& plane, const PxVec3& point)
        {
            return plane.m_origin + plane.m_first_axis * point[0] + plane.m_second_axis * point[1] + plane.m_normal * point[2];
        }

        static inline PxVec3 cartesian_coordinates(const plane& plane, const PxVec2& point)
        {
            return plane.m_origin + plane.m_first_axis * point[0] + plane.m_second_axis * point[1];
        }

    public:
        friend bool operator==(const plane& lhs, const plane& rhs);
        friend bool operator!=(const plane& lhs, const plane& rhs);

    private:
        PxVec3 m_origin;
        PxVec3 m_first_axis;
        PxVec3 m_second_axis;
        PxVec3 m_normal;
    };

    inline bool operator==(const plane& lhs, const plane& rhs)
    {
        return lhs.m_origin == rhs.m_origin &&
            lhs.m_first_axis == rhs.m_first_axis &&
            lhs.m_second_axis == rhs.m_second_axis &&
            lhs.m_normal == rhs.m_normal;
    }

    inline bool operator!=(const plane& lhs, const plane& rhs)
    {
        return lhs.m_origin != rhs.m_origin ||
            lhs.m_first_axis != rhs.m_first_axis ||
            lhs.m_second_axis != rhs.m_second_axis ||
            lhs.m_normal != rhs.m_normal;
    }

    class curve
    {
    public:
        virtual float projected_length(const PxVec3& point) const = 0;
        virtual void frame(float s, PxTransform& result) const = 0;
    };

    class circle : public curve
    {
    public:
        inline circle(plane plane, float radius)
        : m_plane(plane),
          m_radius(radius)
        {
            // Nothing to do.
        }

    public:
        // The circumference.
        inline float length() const
        {
            return m_radius * PxTwoPi;
        }

        // Evaluates the frame on the circle at the specified arc-length.
        inline void frame(float s, PxTransform& result) const
        {
            PxVec3 position = this->point(s);
            PxVec3 tangent = this->tangent(s);
            PxVec3 normal = this->normal(s);
            PxVec3 bitangent = tangent.cross(normal);

            PxQuat& q = result.q;
            result.p = position;

            float m0 = tangent.x;
            float m1 = tangent.y;
            float m2 = tangent.z;
            float m4 = normal.x;
            float m5 = normal.y;
            float m6 = normal.z;
            float m8 = bitangent.x;
            float m9 = bitangent.y;
            float m10 = bitangent.z;

            // Remove the scale from the matrix
            float lx = m0 * m0 + m1 * m1 + m2 * m2;
            float ly = m4 * m4 + m5 * m5 + m6 * m6;
            float lz = m8 * m8 + m9 * m9 + m10 * m10;
            if (lx <= 0.0F || ly <= 0.0F || lz <= 0.0F) { result = PxTransform(PxIdentity); }
            lx = 1.0F / sqrtf(lx);
            ly = 1.0F / sqrtf(ly);
            lz = 1.0F / sqrtf(lz);
            m0 *= lx; m1 *= lx; m2 *= lx;
            m4 *= ly; m5 *= ly; m6 *= ly;
            m8 *= lz; m9 *= lz; m10 *= lz;

            const float trace = m0 + m5 + m10;
            if (trace > 0.0F)
            {
                const float S = sqrtf(trace + 1.0F);
                q.w = S * 0.5F;
                const float t = 0.5F / S;
                q.x = (m6 - m9) * t;
                q.y = (m8 - m2) * t;
                q.z = (m1 - m4) * t;
                return;
            }

            if (m0 >= m5 && m0 >= m10)
            {
                const float S = sqrtf(1.0F + m0 - m5 - m10);
                const float invs = 0.5F / S;
                q.x = 0.5f * S;
                q.y = (m1 + m4) * invs;
                q.z = (m2 + m8) * invs;
                q.w = (m6 - m9) * invs;
                return;
            }
            else if (m5 > m10)
            {
                const float S = sqrtf(1.0F + m5 - m0 - m10);
                const float invs = 0.5F / S;
                q.x = (m4 + m1) * invs;
                q.y = 0.5F * S;
                q.z = (m9 + m6) * invs;
                q.w = (m8 - m2) * invs;
                return;
            }
            else
            {
                const float S = sqrtf(1.0F + m10 - m0 - m5);
                const float invs = 0.5F / S;
                q.x = (m8 + m2) * invs;
                q.y = (m9 + m6) * invs;
                q.z = 0.5F * S;
                q.w = (m1 - m4) * invs;
            }
        }

        // Evaluates the point on the circle at the specified arc-length.
        PxVec3 point(float s) const
        {
            const float r = m_radius;
            if (r <= 0) { return PxVec3(0); }
            const float theta = s / r;
            const float st = PxSin(theta);
            const float ct = PxCos(theta);
            const PxVec3 o = m_plane.origin();
            const PxVec3 n = m_plane.normal();
            const PxVec3 f = m_plane.first_axis();
            const float sax = n.y * f.z - n.z * f.y;
            const float say = n.z * f.x - n.x * f.z;
            const float saz = n.x * f.y - n.y * f.x;
            const float x = o.x + f.x * r * ct + sax * r * st;
            const float y = o.y + f.y * r * ct + say * r * st;
            const float z = o.z + f.z * r * ct + saz * r * st;
            return PxVec3(x, y, z);
        }

        // Evaluates the tangent to the circle at the specified arc-length.
        PxVec3 tangent(float s) const
        {
            const float r = m_radius;
            if (r <= 0) { return PxVec3(0); }
            const PxVec3 n = m_plane.normal();
            const PxVec3 f = m_plane.first_axis();
            const float sax = n.y * f.z - n.z * f.y;
            const float say = n.z * f.x - n.x * f.z;
            const float saz = n.x * f.y - n.y * f.x;
            const float theta = s / r;
            const float st = PxSin(theta);
            const float ct = PxCos(theta);
            const float x = sax * ct - f.x * st;
            const float y = say * ct - f.y * st;
            const float z = saz * ct - f.z * st;
            const float l = sqrtf(x * x + y * y + z * z);
            return PxVec3(x / l, y / l, z / l);
        }

        // Evaluates the normal to the circle at the specified arc-length.
        PxVec3 normal(float s) const
        {
            const float r = m_radius;
            if (r <= 0) { return PxVec3(0); }
            const PxVec3 n = m_plane.normal();
            const PxVec3 f = m_plane.first_axis();
            const float sax = n.y * f.z - n.z * f.y;
            const float say = n.z * f.x - n.x * f.z;
            const float saz = n.x * f.y - n.y * f.x;
            const float theta = s / r;
            const float st = PxSin(theta);
            const float ct = PxCos(theta);
            const float x = -f.x * ct - sax * st;
            const float y = -f.y * ct - say * st;
            const float z = -f.z * ct - saz * st;
            const float l = sqrtf(x * x + y * y + z * z);
            return PxVec3(x / l, y / l, z / l);
        }

        // Calculates the arc-length for the specified point when projected onto the circle
        float projected_length(const PxVec3& point) const
        {
            return m_radius * plane::projected_angle(m_plane, point);
        }

    private:
        plane m_plane;
        float m_radius;
    };

}

class PulleyJoint : public physx::PxConstraintConnector
{
public:
    struct PathJointData : physx::Ext::JointData
	{
        const physx::curve* path;
        physx::PxD6JointDrive drive_settings;
        float drive_position;
        float drive_velocity;
		physx::PxJointLinearLimitPair limit;
		physx::PxPathJointFlags joint_flags;

		PathJointData(const physx::curve* path, const physx::PxJointLinearLimitPair& pair)
        : path(path),
          drive_settings(),
          drive_position(0),
          drive_velocity(0),
          limit(pair),
          joint_flags(0)
        {
            // Nothing to do.
        }
	};

	static const physx::PxU32 TYPE_ID = physx::PxConcreteType::eFIRST_USER_EXTENSION;

	PulleyJoint(physx::PxPhysics& physics, const physx::curve* path, 
				physx::PxRigidBody& body0, const physx::PxTransform& localFrame0,
			    physx::PxRigidBody& body1, const physx::PxTransform& localFrame1);

	void release();

	// attribute accessor and mutators

	//void			setAttachment0(const physx::PxVec3& pos);
	//physx::PxVec3	getAttachment0() const;

	//void			setAttachment1(const physx::PxVec3& pos);
	//physx::PxVec3	getAttachment1() const;

	//void			setDistance(physx::PxReal totalDistance);
	//physx::PxReal	getDistance() const;
	
	//void			setRatio(physx::PxReal ratio);
	//physx::PxReal	getRatio() const;

	// PxConstraintConnector boilerplate

    void            setFlags(physx::PxPathJointFlags flags) { m_data.joint_flags = flags; }
    void            setDrivePosition(float position) { m_data.drive_position = position; }
    void            setDriveVelocity(float velocity) { m_data.drive_velocity = velocity; }
    void            setDriveSettings(const physx::PxD6JointDrive& settings) { m_data.drive_settings = settings; }

	void*			prepareData();
	void			onConstraintRelease();
	void			onComShift(physx::PxU32 actor);
	void			onOriginShift(const physx::PxVec3& shift);
	void*			getExternalReference(physx::PxU32& typeID);

	bool			updatePvdProperties(physx::pvdsdk::PvdDataStream&,
										const physx::PxConstraint*,
										physx::PxPvdUpdateType::Enum) const { return true; }
	void			updateOmniPvdProperties() const { }
	physx::PxBase*	getSerializable() { return NULL; }

	virtual physx::PxConstraintSolverPrep getPrep() const;

	virtual const void* getConstantBlock() const { return &m_data; }

	physx::PxRigidBody*		m_bodies[2];
	physx::PxTransform		m_local_poses[2];

	physx::PxConstraint*	m_constraint;
	PathJointData			m_data;

	~PulleyJoint() {}
};

#endif

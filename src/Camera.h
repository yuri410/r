#pragma once

#include "Common.h"

struct Camera
{
public:
	Camera(Vec3 pos, float fovY, float aspectRatio);

	Ray ComputeCameraRay(float x, float y) const;
private:
	float m_fovY;

	float m_aspectRatio;

	GL::Mat4 m_viewProj;
	GL::Mat4 m_invViewProj;
};
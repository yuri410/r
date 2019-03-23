#include "Camera.h"
#include "Math.h"

Vec3 Unproject(const Vec3& p, const GL::Mat4& invViewProj)
{
	Vec3 v = p * 2 - Vec3(1.0f, 1.0f, 1.0f);

	Vec3 res = invViewProj * v;
	float w = invViewProj.m[3] * v.X + invViewProj.m[7] * v.Y + invViewProj.m[11] * v.Z + invViewProj.m[15];
	return res / w;
}

Camera::Camera(Vec3 pos, float fovY, float aspectRatio)
	: m_fovY(fovY), m_aspectRatio(aspectRatio)
{
	GL::Mat4 view = GL::Mat4::LookAt(pos, GL::Vec3(0, 0, 0), GL::Vec3(0, 0, 1));
	GL::Mat4 proj = GL::Mat4::Perspective(GL::Rad(fovY), aspectRatio, 0.1f, 1000.0f);

	m_viewProj = proj*view;
	m_invViewProj = m_viewProj.Inverse();
}

Ray Camera::ComputeCameraRay(float x, float y) const
{
	Vec3 startPt = Unproject(Vec3(x, y, 0), m_invViewProj);
	Vec3 endPt = Unproject(Vec3(x, y, 1), m_invViewProj);

	Ray res = { startPt, (endPt-startPt).Normal() };
	return res;
}
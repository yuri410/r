#pragma once

#include "Common.h"

struct SceneIntersection
{
	Vec3 m_point;
	Vec3 m_normal;
	SceneObject* m_sceneObject = nullptr;
};

class Scene
{
public:
	Scene(float cameraAspectRatio, float cameraFovY);
	~Scene();

	bool Intersect(const Ray& ray, SceneIntersection& intersection);

	Camera* GetCamera() const { return m_camera; }

private:

	Camera* m_camera = nullptr;
	
	std::vector<SceneObject*> m_sceneObjects;
};
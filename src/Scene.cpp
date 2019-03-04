#include "Scene.h"
#include "Camera.h"
#include "Geometry.h"

Scene::Scene(float cameraAspectRatio, float cameraFovY)
{
	m_camera = new Camera(Vec3(3, 3, 3), cameraFovY, cameraAspectRatio);


}

Scene::~Scene()
{
	for (auto& obj : m_sceneObjects)
	{
		delete obj; obj = nullptr;
	}
	m_sceneObjects.clear();
}

bool Scene::Intersect(const Ray& ray, SceneIntersection& intersection)
{
	//SphereShape sphere(Vec3(0, 0, 1), 1);
	//
	//Vec3 pt, n;
	//return sphere.Intersects(ray, pt, n);
	
	QuadShape quad(Vec3(0, 0, 0), Vec3(0, 0, 1), Vec3(0, 1, 0), 3, 1);
	Vec3 pt, n;
	return quad.Intersects(ray, pt, n);
}
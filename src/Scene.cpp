#include "Scene.h"
#include "Camera.h"
#include "Geometry.h"
#include "SceneObject.h"
#include "Material.h"

Scene::Scene(float cameraAspectRatio, float cameraFovY)
{
	m_camera = new Camera({ 13, 2, 1 }, cameraFovY, cameraAspectRatio);

	float rs = 8;

	//m_sceneObjects.push_back(new SceneObject(new Material({ 0.5f,0.5f,0.5f }, { 100,100,100 }, 0), new SphereShape({ 0, 0,  4.8f }, 1)));
	m_sceneObjects.push_back(new SceneObject(new Material(Vec3{ 0.5f,0.5f,0.5f }), new SphereShape({ -1, -2.5f, -1.5f }, 1.5f)));
	m_sceneObjects.push_back(new SceneObject(new Material(Vec3{ 0.5f,0.5f,0.5f }), new SphereShape({ -2, 1, -2.5f }, 1.5f)));

	m_sceneObjects.push_back(new SceneObject(new Material(Vec3{ 0.75f,0.75f,0.75f }), new QuadShape({ -rs*0.5f, 0, 0 }, {  1, 0,  0 }, { 0,  1, 0 }, rs, rs)));
	//m_sceneObjects.push_back(new SceneObject(new Material(Vec3{ 0.75f,0.75f,0.75f }), new QuadShape({  rs*0.5f, 0, 0 }, { -1, 0,  0 }, { 0, -1, 0 }, rs, rs)));
	m_sceneObjects.push_back(new SceneObject(new Material(Vec3{ 0.75f,0.25f,0.25f }), new QuadShape({ 0,-rs*0.5f, 0 }, { 0,  1, 0 }, { 0, 0,  1 }, rs, rs)));
	m_sceneObjects.push_back(new SceneObject(new Material(Vec3{ 0.25f,0.25f,0.75f }), new QuadShape({ 0, rs*0.5f, 0 }, { 0, -1, 0 }, { 0, 0, -1 }, rs, rs)));
	m_sceneObjects.push_back(new SceneObject(new Material(Vec3{ 0.75f,0.75f,0.75f }), new QuadShape({ 0, 0,-rs*0.5f }, { 0, 0,  1 }, {  1, 0, 0 }, rs, rs)));
	m_sceneObjects.push_back(new SceneObject(new Material(Vec3{ 0.75f,0.75f,0.75f }, Vec3(10,10,10), 0), new QuadShape({ 0, 0, rs*0.5f }, { 0, 0, -1 }, { -1, 0, 0 }, rs, rs)));
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
	
	float minDist = FLT_MAX;
	bool  intersected = false;
	Vec3  pt, n;
	for (SceneObject* so : m_sceneObjects)
	{
		IGeometry* geom = so->GetGeometry();

		if (geom && geom->Intersects(ray, pt, n))
		{
			float distSq = (pt - ray.m_origin).LengthSqr();
			if (distSq < minDist)
			{
				intersected = true;
				minDist = distSq;

				intersection.m_point = pt;
				intersection.m_normal = n;
				intersection.m_sceneObject = so;
			}
		}
	}
	return intersected;
}
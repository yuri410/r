#include "Scene.h"
#include "Camera.h"
#include "Geometry.h"
#include "SceneObject.h"
#include "Material.h"
#include "Math.h"

Scene::Scene(float cameraAspectRatio, float cameraFovY)
{
	const int sceneId = 1;

	if (sceneId == 0)
	{
		m_camera = new Camera({ 13, 2, 1 }, cameraFovY, cameraAspectRatio);

		float rs = 8;

		//m_sceneObjects.push_back(new SceneObject(new Material({ 0.5f,0.5f,0.5f }, { 100,100,100 }, 0), new SphereShape({ 0, 0,  4.8f }, 1)));
		m_sceneObjects.push_back(new SceneObject(new Material(Vec3{ 0.5f,0.5f,0.5f }, {}, 1, MaterialType::Diffuse), 
			new CylinderShape({ 0, 0, -2.0f }, { -1,1,1 }, 4, 1)));
		//m_sceneObjects.push_back(new SceneObject(new Material(Vec3{ 0.5f,0.5f,0.5f }), new SphereShape({ -2, -2.5f, -1.5f }, 1.5f)));
		//m_sceneObjects.push_back(new SceneObject(new Material(Vec3{ 0.5f,0.5f,0.5f }, {}, 1, MaterialType::MicrofacetBRDF), new SphereShape({ -2, 1, -2.5f }, 1.5f)));

		//m_sceneObjects.push_back(new SceneObject(new Material(Vec3{ 0.5f,0.5f,0.5f }, {30,30,30},0), new CylinderShape({ 0, 0, rs*0.5f - 0.01f }, { 0,0,1 }, 0.01f, 2)));


		m_sceneObjects.push_back(new SceneObject(new Material(Vec3{ 0.75f,0.75f,0.75f }), new QuadShape({ -rs*0.5f, 0, 0 }, { 1, 0,  0 }, { 0,  1, 0 }, rs, rs)));
		//m_sceneObjects.push_back(new SceneObject(new Material(Vec3{ 0.75f,0.75f,0.75f }), new QuadShape({  rs*0.5f, 0, 0 }, { -1, 0,  0 }, { 0, -1, 0 }, rs, rs)));
		m_sceneObjects.push_back(new SceneObject(new Material(Vec3{ 0.75f,0.25f,0.25f }), new QuadShape({ 0,-rs*0.5f, 0 }, { 0,  1, 0 }, { 0, 0,  1 }, rs, rs)));
		m_sceneObjects.push_back(new SceneObject(new Material(Vec3{ 0.25f,0.25f,0.75f }), new QuadShape({ 0, rs*0.5f, 0 }, { 0, -1, 0 }, { 0, 0, -1 }, rs, rs)));
		m_sceneObjects.push_back(new SceneObject(new Material(Vec3{ 0.75f,0.75f,0.75f }), new QuadShape({ 0, 0,-rs*0.5f }, { 0, 0,  1 }, { 1, 0, 0 }, rs, rs)));
		m_sceneObjects.push_back(new SceneObject(new Material(Vec3{ 0.75f,0.75f,0.75f }, Vec3(10,10,10), 1), new QuadShape({ 0, 0, rs*0.5f }, { 0, 0, -1 }, { -1, 0, 0 }, rs, rs)));
	}
	else if (sceneId == 1)
	{
		int xcount = 4;
		int ycount = 3;

		float xSpan = 24;
		float ySpan = xSpan / xcount*ycount;
		
		cameraFovY = 30;

		float viewDistance = ySpan*0.5f / tan(GL::Rad(cameraFovY / 2));
		m_camera = new Camera({ viewDistance, 0, 0 }, cameraFovY, cameraAspectRatio);

		float cellSize = xSpan / xcount;
		float ballRadius = cellSize*0.5f*0.67f;
		float lightRadius = cellSize*0.5f*0.1f;

		int counter = 0;
		int totalCount = xcount*ycount;
		for (int y = 0; y < ycount; y++)
		{
			for (int x = 0; x < xcount; x++)
			{
				float cellX = (x+0.5f) * cellSize - xSpan / 2;
				float cellY = (ycount - y - 1 + 0.5f) * cellSize - ySpan / 2;
				//float cellY = (y + 0.5f) * cellSize - ySpan / 2;

				Vec3 cellPos = { -cellSize*0.5f, cellX, cellY };

				m_sceneObjects.push_back(new SceneObject(
					new Material(Vec3{ 0.5f,0.5f,0.5f }, {}, max(1e-3,(float)counter / (totalCount - 1)), MaterialType::MicrofacetBRDF),
					new SphereShape(cellPos, ballRadius)
				));

				float ws = cellSize - 0.02f;
				m_sceneObjects.push_back(new SceneObject(new Material(Vec3{ 0.75f,0.75f,0.75f }), new QuadShape(Vec3{ -ws*0.5f, 0, 0 }+cellPos, { 1, 0,  0 }, { 0,  1, 0 }, ws, ws)));
				m_sceneObjects.push_back(new SceneObject(new Material(Vec3{ 0.75f,0.50f,0.50f }), new QuadShape(Vec3{ 0,-ws*0.5f, 0 }+cellPos, { 0,  1, 0 }, { 0, 0,  1 }, ws, ws)));
				m_sceneObjects.push_back(new SceneObject(new Material(Vec3{ 0.50f,0.75f,0.75f }), new QuadShape(Vec3{ 0, ws*0.5f, 0 }+cellPos, { 0, -1, 0 }, { 0, 0, -1 }, ws, ws)));
				m_sceneObjects.push_back(new SceneObject(new Material(Vec3{ 0.50f,0.75f,0.50f }), new QuadShape(Vec3{ 0, 0,-ws*0.5f }+cellPos, { 0, 0,  1 }, { 1, 0, 0 }, ws, ws)));
				m_sceneObjects.push_back(new SceneObject(new Material(Vec3{ 0.75f,0.75f,0.50f }, {5,5,5}, 1), new QuadShape(Vec3{ 0, 0, ws*0.5f }+cellPos, { 0, 0, -1 }, { -1, 0, 0 }, ws, ws)));

				
				{
					m_sceneObjects.push_back(new SceneObject(
						new Material(Vec3{ 0.5f,0.5f,0.5f }, { 150,150,150 }, 0),
						new DiskShape(Vec3{ 0,-cellSize*0.5f,-ws*0.5f + 0.01f }+cellPos, { 0,0,1 }, lightRadius))
					);
				}
				if (x == xcount - 1)
				{
					m_sceneObjects.push_back(new SceneObject(
						new Material(Vec3{ 0.5f,0.5f,0.5f }, { 150,150,150 }, 0),
						new DiskShape(Vec3{ 0,cellSize*0.5f,-ws*0.5f + 0.01f }+cellPos, { 0,0,1 }, lightRadius))
					);
				}
				counter++;
			}
		}
	}
	

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
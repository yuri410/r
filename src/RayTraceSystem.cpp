#include "RayTraceSystem.h"
#include "Math.h"
#include "Scene.h"
#include "SceneObject.h"
#include "Material.h"
#include "Camera.h"
#include "BRDF.h"

RayTraceSystem::RayTraceSystem(int imageWidth, int imageHeight, int tileSize, int numThreads)
	: m_imageWidth(imageWidth), m_imageHeight(imageHeight), m_tileSize(tileSize), m_threadCount(numThreads)
{
	m_resultBuffer = new Vec3[imageHeight*imageWidth]();
	m_resultBufferTemp = new Vec3[imageHeight*imageWidth]();
	m_sampleCounter = new int[imageHeight*imageWidth]();
	m_presentTexture.SetFilters(GL::Filter::Nearest, GL::Filter::Nearest);

	m_scene = new Scene(imageWidth / (float)imageHeight, 60);
}

RayTraceSystem::~RayTraceSystem()
{
	for (auto& thread : m_renderThreads)
	{
		delete thread;
	}
	m_renderThreads.clear();

	delete m_scene;
	m_scene = nullptr;
	delete[] m_sampleCounter;
	delete[] m_resultBuffer;
	delete[] m_resultBufferTemp;
}

void RayTraceSystem::Start()
{
	for (int i = 0; i < m_threadCount; i++)
	{
		m_renderThreads.push_back(new RenderThread(this));
	}

	Random rnd;
	bool tiled = true;
	
	if (!tiled)
	{
		std::vector<int> pixelIndices;

		for (int i = 0; i < m_imageHeight; i++)
		{
			for (int j = 0; j < m_imageWidth; j++)
			{
				pixelIndices.push_back(i*m_imageWidth + j);
			}
		}
		for (size_t i = 0; i < pixelIndices.size(); i++)
		{
			int index = rnd.IntSample() % (i+1);

			std::swap(pixelIndices[index], pixelIndices[i]);
		}

		int tid = 0;
		for (int pidx : pixelIndices)
		{
			int x = pidx % m_imageWidth;
			int y = pidx / m_imageWidth;

			RenderTask task = { x,y };
			m_renderThreads[tid++]->AddJob(task);
			tid = tid % m_threadCount;
		}
	}
	else
	{
		int tileSize = m_tileSize;

		for (int i = 0; i < m_imageHeight; i += tileSize)
		{
			for (int j = 0; j < m_imageWidth; j += tileSize)
			{
				int maxX = std::min<int>(m_imageWidth, j + tileSize);
				int maxY = std::min<int>(m_imageHeight, i + tileSize);

				for (int m = i; m < maxY; m++)
				{
					for (int n = j; n < maxX; n++)
					{
						RenderTask task = { n,m };
						int threadId = m_threadCount > 1 ? (rnd.IntSample() % m_threadCount) : 0;
						m_renderThreads[threadId]->AddJob(task);
					}
				}
			}
		}
	}

}

void RayTraceSystem::CopyBufferData()
{
	m_resultMutex.lock();
	memcpy(m_resultBufferTemp, m_resultBuffer, sizeof(Vec3)*m_imageHeight*m_imageWidth);
	m_resultMutex.unlock();

	m_presentTexture.Image2D(m_resultBufferTemp, GL::DataType::Float, GL::Format::RGB, m_imageWidth, m_imageHeight, GL::InternalFormat::RGB32F);
}

void RayTraceSystem::ReportResult(int x, int y, Vec3 p)
{
	if (isnan(p.X) || isnan(p.Y) || isnan(p.Z))
		return;

	m_resultMutex.lock();
	
	int idx = y*m_imageWidth + x;
	int sampleCount = m_sampleCounter[idx];
	m_sampleCounter[idx]++;

	if (sampleCount == 0)
	{
		m_resultBuffer[idx] = p;
	}
	else
	{
		m_resultBuffer[idx] = m_resultBuffer[idx] * ((float)sampleCount / (sampleCount + 1)) + p / (sampleCount + 1);
	}

	m_resultMutex.unlock();
}


Vec3 RayTraceSystem::Trace(const Ray& ray, int depth, Random& rnd)
{
	const float bias = 0.00001f;
	//const float bias = 0.01f;
	
	Vec3 result;
	SceneIntersection intersection;
	if (m_scene->Intersect(ray, intersection))
	{
		Material* mtrl = intersection.m_sceneObject->GetMaterial();
		result += mtrl->m_emissive;

		if (depth < 5)
		{
			if (mtrl->m_type == MaterialType::Diffuse)
			{
				float e1 = rnd.FloatSample();
				float e2 = rnd.FloatSample();

				float pdf;
				HemisphereCoord c;
				ImportanceSampleCosine(e1, e2, c, pdf);

				Vec3 xi = c.GetCartesianCoord_AutoFrameZ(intersection.m_normal);

				float ndl = c.m_cosTheta;

				Ray nextRay = { intersection.m_point + intersection.m_normal * bias, xi };
				result += Trace(nextRay, depth + 1, rnd) * mtrl->m_albedo * (ndl / (PI*pdf));
			}
			else if (mtrl->m_type == MaterialType::Reflection)
			{
				float ndr = intersection.m_normal.Dot(ray.m_direction);

				Vec3 r = (ray.m_direction - intersection.m_normal * (2 * ndr)).Normal();

				Ray nextRay = { intersection.m_point + intersection.m_normal * bias, r };
				result += Trace(nextRay, depth + 1, rnd) * mtrl->m_reflectivity;
			}
			else if (mtrl->m_type == MaterialType::MicrofacetBRDF)
			{
				float e1 = rnd.FloatSample();
				float e2 = rnd.FloatSample();

				float roughness = mtrl->m_roughness;

				float pdf;
				HemisphereCoord c;
				ImportanceSampleGGX(e1, e2, roughness, c, pdf);

				Vec3 h = c.GetCartesianCoord_AutoFrameZ(intersection.m_normal);

				float hdr = h.Dot(ray.m_direction);

				Vec3 xi = (ray.m_direction - h * (2 * hdr)).Normal();
				
				float ndl = xi.Dot(intersection.m_normal);
				float ndh = (h.Dot(intersection.m_normal)); // cosTheta;
				float ndv = -(ray.m_direction.Dot(intersection.m_normal));
				float ldh = (h.Dot(xi));
				float hdv = abs(h.Dot(ray.m_direction));

				pdf = pdf / (4 * hdv);

				if (ndv > 0 && ndl > 0)
				{
					float D = BRDF::NDF(ndh, roughness);
					float G = BRDF::GSF(ndl, ndv, roughness);
					Vec3 Fr = BRDF::F(mtrl->m_reflectivity, ldh);

					Ray nextRay = { intersection.m_point + intersection.m_normal * bias, xi };
					result += Trace(nextRay, depth + 1, rnd) * BRDF::CookTorrance(Fr, D, G, ndl, ndv) * (ndl / pdf);
				}
			}
		}
	}
	return result;
}

//////////////////////////////////////////////////////////////////////////

RayTraceSystem::RenderThread::RenderThread(RayTraceSystem* sys)
	: m_thread(ThreadEntry, this), m_system(sys), m_terminating(false), m_isIdle(false)
{

}

RayTraceSystem::RenderThread::~RenderThread()
{
	m_queueMutex.lock();
	while (!m_jobQueue.empty())
		m_jobQueue.pop();
	m_queueMutex.unlock();

	m_terminating = true;
	m_queueWaiter.notify_all();
	if (m_thread.joinable())
		m_thread.join();
}

void RayTraceSystem::RenderThread::Run()
{
	Random rnd;

	float invWidth = 1.0f / m_system->m_imageWidth;
	float invHeight = 1.0f / m_system->m_imageHeight;

	while (!m_terminating)
	{
		volatile bool hasTask = false;
		RenderTask task;
		{
			std::unique_lock<std::mutex> lock(m_queueMutex);
			if (m_jobQueue.size() == 0)
			{
				m_isIdle = true;

				m_queueWaiter.wait(lock);
				m_isIdle = false;
			}
			else
			{
				task = m_jobQueue.front();
				m_jobQueue.pop();
				hasTask = true;
			}
		}

		if (m_terminating)
			break;

		if (hasTask)
		{
			assert(task.m_x < m_system->m_imageWidth);
			assert(task.m_y < m_system->m_imageHeight);

			const int sampleCount = 16;

			Vec3 result;
			for (int s = 0; s < sampleCount; s++)
			{
				float r1 = rnd.FloatSample() * 2 - 1;
				float r2 = rnd.FloatSample() * 2 - 1;

				float dx = (1 - sqrt(fabs(r1))) * (r1 > 0 ? 1 : -1) + 0.5f;
				float dy = (1 - sqrt(fabs(r2))) * (r2 > 0 ? 1 : -1) + 0.5f;

				float x = (task.m_x + dx) * invWidth;
				float y = (task.m_y + dy) * invHeight;

				if (m_system->m_scene && m_system->m_scene->GetCamera())
				{
					Ray ray = m_system->m_scene->GetCamera()->ComputeCameraRay(x, y);
			
					result += m_system->Trace(ray, 0, rnd);
				}
			}

			m_system->ReportResult(task.m_x, task.m_y, result / sampleCount);
			AddJob(task);
		}
	}
}

void RayTraceSystem::RenderThread::ThreadEntry(void* renderThread)
{
	if (renderThread)
	{
		RenderThread* obj = reinterpret_cast<RenderThread*>(renderThread);
		obj->Run();
	}
}

void RayTraceSystem::RenderThread::AddJob(const RenderTask& task)
{
	m_queueMutex.lock();
	m_jobQueue.push(task);

	m_queueWaiter.notify_all();
	m_queueMutex.unlock();
}
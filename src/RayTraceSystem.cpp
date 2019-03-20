#include "RayTraceSystem.h"
#include "Random.h"
#include "Scene.h"
#include "SceneObject.h"
#include "Material.h"
#include "Camera.h"

RayTraceSystem::RayTraceSystem(int imageWidth, int imageHeight, int tileSize, int numThreads)
	: m_imageWidth(imageWidth), m_imageHeight(imageHeight), m_tileSize(tileSize), m_threadCount(numThreads)
{
	m_resultBuffer = new Vec3[imageHeight*imageWidth]();
	m_resultBufferTemp = new Vec3[imageHeight*imageWidth]();
	m_presentTexture.SetFilters(GL::Filter::Nearest, GL::Filter::Nearest);

	m_samples = new ComposingPixel[imageHeight*imageWidth]();

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
	ComposingPixel& cp = m_samples[y*m_imageWidth + x];
	cp.m_samples.push_back(p);

	Vec3 res;
	for (auto& s : cp.m_samples)
	{
		res += s;
	}
	if (cp.m_samples.size())
		res = res / (float)cp.m_samples.size();

	m_resultMutex.lock();
	m_resultBuffer[y*m_imageWidth + x] = res;
	m_resultMutex.unlock();
}

float sqr(float x) { return x*x; }

// GGX
float NDF(float NdotH, float roughness)
{
	float pi = 3.1415927f;

	float a2 = roughness*roughness;
	return a2 / (pi * sqr((a2 - 1)*NdotH*NdotH + 1));
}

// Schlick 
float GSF(float NdotL, float NdotV, float roughness)
{
	float k = roughness / 2;

	float SmithL = (NdotL) / (NdotL * (1 - k) + k);
	float SmithV = (NdotV) / (NdotV * (1 - k) + k);

	return (SmithL * SmithV);
}

Vec3 F(Vec3 f0, float LdotH)
{
	return f0 + (Vec3(1.0f, 1.0f, 1.0f) - f0) * powf(LdotH, 5);
}


Vec3 RayTraceSystem::Trace(const Ray& ray, int depth, Random& rnd)
{
	const float pi = 3.14159f;
	const float bias = 0.01f;

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
				Vec3 basisX, basisY, basisZ = intersection.m_normal;
				muAutoFrameZ(basisZ, basisX, basisY);
				
				bool uniform = false;
				float theta = rnd.FloatSample() * 2 * pi;
				float sinPhi, cosPhi;
				float pdf;

				if (uniform)
				{
					float r2 = rnd.FloatSample();
					cosPhi = r2;
					sinPhi = sqrtf(1 - cosPhi*cosPhi);
					pdf = 1.0f / (2 * pi);
				}
				else
				{
					float r2 = rnd.FloatSample();
					sinPhi = sqrtf(r2);
					cosPhi = sqrtf(1 - r2);
					pdf = cosPhi / pi;
				}

				Vec3 xi = basisX*(cos(theta)*sinPhi) + basisY*(sin(theta)*sinPhi) + basisZ * cosPhi;
				xi = xi.Normal();

				float ndl = cosPhi;

				Ray ray2 = { intersection.m_point + basisZ * bias, xi };
				result += Trace(ray2, depth + 1, rnd) * mtrl->m_albedo * (ndl / (pi*pdf));
			}
			else if (mtrl->m_type == MaterialType::Reflection)
			{
				float ndr = intersection.m_normal.Dot(ray.m_direction);

				Vec3 r = ray.m_direction - intersection.m_normal * (2 * ndr);
				r = r.Normal();

				Ray ray2 = { intersection.m_point + intersection.m_normal * bias, r };
				result += Trace(ray2, depth + 1, rnd);// * mtrl->m_albedo;
			}
			else if (mtrl->m_type == MaterialType::MicrofacetBRDF)
			{
				Vec3 basisX, basisY, basisZ = intersection.m_normal;
				muAutoFrameZ(basisZ, basisX, basisY);

				int   mode = 2;
				float theta = rnd.FloatSample() * 2 * pi;
				float sinPhi, cosPhi;
				float pdf;

				float roughness = mtrl->m_roughness;

				if (mode == 0)
				{
					// uniform
					float rn = rnd.FloatSample();
					cosPhi = rn;
					sinPhi = sqrtf(1 - cosPhi*cosPhi);
					pdf = 1.0f / (2 * pi);
				}
				else if (mode == 1)
				{
					// cosine
					float rn = rnd.FloatSample();
					sinPhi = sqrtf(rn);
					cosPhi = sqrtf(1 - rn);
					pdf = cosPhi / pi;
				}
				else if (mode == 2)
				{
					// ggx
					float rn = rnd.FloatSample();

					float a2 = roughness*roughness;

					float rc = (1 - rn) / ((a2 - 1.0) *rn + 1);
					rc = min(1.0f, max(0.0f, rc));

					cosPhi = sqrtf(rc);
					sinPhi = sqrtf(1 - rc);

					pdf = a2 * cosPhi / (pi * sqr((a2 - 1)*cosPhi*cosPhi + 1));
				}

				Vec3 xi, h;
				{
					h = basisX*(cos(theta)*sinPhi) + basisY*(sin(theta)*sinPhi) + basisZ * cosPhi;
					h = h.Normal();

					float hdr = h.Dot(ray.m_direction);

					xi = ray.m_direction - h * (2 * hdr);
					xi.Normal();
				}

				float ndl = xi.Dot(intersection.m_normal);
				float ndh = (h.Dot(intersection.m_normal)); // cosPhi;
				float ndv = -(ray.m_direction.Dot(intersection.m_normal));
				float ldh = (h.Dot(xi));
				float hdv = abs(h.Dot(ray.m_direction));

				pdf = pdf / (4 * hdv);

				if (ndv>0 && ndl>0 )
				{
					float D = NDF(ndh, roughness);
					float G = GSF(ndl, ndv, roughness);
					Vec3 Fr = F(mtrl->m_reflectivity, ldh);

					Ray ray2 = { intersection.m_point + intersection.m_normal * bias, xi };
					result += Trace(ray2, depth+1, rnd) * Fr * (D * G) / (4 * ndl * ndv) * (ndl / pdf);
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

			for (int s = 0; s < 64; s++)
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
			
					Vec3 result = m_system->Trace(ray, 0, rnd);
					m_system->ReportResult(task.m_x, task.m_y, result);
				}

			}
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
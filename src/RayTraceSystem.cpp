#include "RayTraceSystem.h"
#include "Random.h"
#include "Scene.h"
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

	int tileSize = m_tileSize;
	Random rnd;
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

			for (int s = 0; s < 4; s++)
			{
				float r1 = rnd.FloatSample() * 2 - 1;
				float r2 = rnd.FloatSample() * 2 - 1;

				float dx = sqrt(fabs(r1)) * (r1 > 0 ? 1 : -1) + 0.5f;
				float dy = sqrt(fabs(r2)) * (r2 > 0 ? 1 : -1) + 0.5f;

				//float r1 = rnd.FloatSample() * 2;
				//float r2 = rnd.FloatSample() * 2;

				//float dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
				//float dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
				//
				//dx += 0.5f; dy += 0.5f;

				float x = (task.m_x + dx) * invWidth;
				float y = (task.m_y + dy) * invHeight;

				if (m_system->m_scene && m_system->m_scene->GetCamera())
				{
					Ray ray = m_system->m_scene->GetCamera()->ComputeCameraRay(x, y);
			
					SceneIntersection intersection;
					if (m_system->m_scene->Intersect(ray, intersection))
					{
						m_system->ReportResult(task.m_x, task.m_y, Vec3(1, 1, 1));
					}
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
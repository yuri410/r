#pragma once

#include "Common.h"
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <atomic>

class RayTraceSystem
{
public:
	RayTraceSystem(int imageWidth, int imageHeight, int tileSize, int numThreads);
	~RayTraceSystem();

	void Start();

	void CopyBufferData();

	const GL::Texture& GetResultTexture() { return m_presentTexture; }
private:
	struct ComposingPixel
	{
		std::vector<Vec3> m_samples;
	};

	struct RenderTask
	{
		int m_x;
		int m_y;
	};

	struct RenderThread
	{
		RayTraceSystem* m_system;

		std::atomic<bool> m_terminating;
		std::atomic<bool> m_isIdle;
		std::mutex m_queueMutex;
		std::condition_variable m_queueWaiter;

		std::queue<RenderTask> m_jobQueue;

		std::thread m_thread;

		RenderThread(RayTraceSystem* sys);
		~RenderThread();

		RenderThread(const RenderThread&) = delete;
		RenderThread(RenderThread&&) = delete;
		RenderThread& operator=(const RenderThread&) = delete;
		RenderThread& operator=(RenderThread&&) = delete;

		static void ThreadEntry(void* rayTraceSys);

		void Run();
		void AddJob(const RenderTask& task);
	};

	void ReportResult(int x, int y, Vec3 p);

	Vec3* m_resultBuffer;
	Vec3* m_resultBufferTemp;
	std::mutex m_resultMutex;

	ComposingPixel* m_samples;

	std::vector<RenderThread*> m_renderThreads;

	int m_imageWidth;
	int m_imageHeight;
	int m_threadCount;
	int m_tileSize;
	GL::Texture m_presentTexture;

	Scene* m_scene = nullptr;
};
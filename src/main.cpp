#include <cmath>
#include "Common.h"
#include "RayTraceSystem.h"

const int WindowWidth = 1280;
const int WindowHeight = 720;

int main()
{
	GL::Window window(WindowWidth, WindowHeight, "Ray Tracer", GL::WindowStyle::Close);
	GL::Context& gl = window.GetContext();

	GL::Shader vert(GL::ShaderType::Vertex, GLSL(
		in vec2 position;
		out vec2 texCoord;
		void main()
		{
			texCoord = (position-vec2(1))*0.5;
			gl_Position = vec4(position, 0.0, 1.0); 
		})
	);

	GL::Shader frag(GL::ShaderType::Fragment, GLSL(
		uniform sampler2D tex;
		uniform float exposure;
		in vec2 texCoord;
		out vec4 color;
		void main()
		{
			vec3 hdrColor = texture(tex, texCoord).rgb;
			color.rgb = vec3(1.0) - exp(-hdrColor * exposure);
			color.a = 1;
		})
	);

	GL::Program program(vert, frag);

	const float vertices[] = {
		-3, -1,
		 1,  3,
		 1, -1
	};
	GL::VertexBuffer vbo(vertices, sizeof(vertices), GL::BufferUsage::StaticDraw);

	GL::VertexArray vao;
	vao.BindAttribute(program.GetAttribute("position"), vbo, GL::Type::Float, 2, 0, 0);
	
	program.SetUniform(program.GetUniform("tex"), 0);
	program.SetUniform(program.GetUniform("exposure"), 1.0f);

	RayTraceSystem rayTracer(WindowWidth, WindowHeight, 16, 8);
	rayTracer.Start();

	GL::Event ev;
	while (window.IsOpen())
	{
		while (window.GetEvent(ev));

		rayTracer.CopyBufferData();

		gl.Clear();

		gl.BindTexture(rayTracer.GetResultTexture(), 0);
		gl.DrawArrays(vao, GL::Primitive::Triangles, 0, 3);

		window.Present();
	}

	return 0;
}

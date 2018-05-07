#version 400

uniform mat4 eye_projection;
uniform mat4 world_to_eye;

layout(location = 0) in vec3 vertex_position;
layout(location = 1) in vec2 vertex_texture;
layout(location = 2) in vec3 vertex_normal;

out vec2 v_texture;
out vec3 v_normal;

void main()
{
    v_texture = vertex_texture;
	v_normal = vertex_normal;
    gl_Position =  eye_projection * world_to_eye * vec4(vertex_position,1.0);
}
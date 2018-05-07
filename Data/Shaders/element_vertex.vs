#version 400

uniform mat4 eye_projection;
uniform mat4 world_to_eye;

layout(location = 0) in vec3 vertex_position;
layout(location = 1) in vec2 vertex_texture;

out vec2 v_texture;

void main()
{
    v_texture = vertex_texture;
    gl_Position =  eye_projection * world_to_eye * vec4(vertex_position,1.0);
}

#version 400

uniform mat4 eye_projection;
uniform mat4 world_to_eye;

layout(location = 0) in vec3 vertex_position;
layout(location = 1) in vec3 vertex_normal;

out vec3 v_eye_to_position;
out vec3 v_normal;

void main()
{
    v_normal = vertex_normal;
    // Camera Clipping Coordinates
	v_eye_to_position = vec3( world_to_eye * vec4(vertex_position,1.0));
    gl_Position =  eye_projection * world_to_eye * vec4(vertex_position,1.0);
}

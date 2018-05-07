#version 400

uniform sampler2D color_texture;
in vec2 v_texture;
layout (location = 0) out vec4 FragColor;

void main(){
   vec3 color_tex = vec3(texture(color_texture, v_texture));
   FragColor = vec4(color_tex,1.0);
}

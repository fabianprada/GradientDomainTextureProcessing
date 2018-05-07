#version 400

uniform vec3 light_direction;

uniform sampler2D color_texture;
in vec2 v_texture;
in vec3 v_normal;
layout (location = 0) out vec4 FragColor;

void main(){
   vec3 surface_normal = normalize(v_normal);
   vec3 position_to_light = normalize(-light_direction);
   float attenuation;

   if(dot(position_to_light,surface_normal) > 0.f){
		//Diffuse attenuation
		attenuation = max(dot(surface_normal,position_to_light), 0.0)*1.5;


		//Specular attenuation
		//vec3 position_to_eye = normalize(-v_eye_to_position);
        //vec3 specular_direction = 2.f*surface_normal*dot(position_to_light,surface_normal) - position_to_light;
        //float specular_attenuation = pow(max(dot(specular_direction,position_to_eye), 0.0),specular_falloff);
		//lighting_color += specular_attenuation*light_specular;
    }

   vec3 color_tex = vec3(texture(color_texture, v_texture))*attenuation;
   FragColor = vec4(color_tex,1.0);
}
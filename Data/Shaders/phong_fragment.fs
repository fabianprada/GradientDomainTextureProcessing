#version 400

uniform vec3 light_direction;
uniform vec3 light_diffuse;
uniform vec3 light_specular;
uniform vec3 light_ambient;
uniform float specular_falloff;

in vec3 v_eye_to_position;
in vec3 v_normal;

layout (location = 0) out vec4 FragColor;

float DiffuseAttenuation()
{
    vec3 position_to_light = normalize(-light_direction);
    vec3 view_normal = normalize(v_normal);
    float attenuation = max(dot(view_normal,position_to_light), 0.0);
    return attenuation;
}

float SpecularAttenuation()
{
    vec3 position_to_eye = normalize(-v_eye_to_position);
    vec3 position_to_light = normalize(-light_direction);
    vec3 view_normal = normalize(v_normal);

    if(dot(position_to_light,view_normal) >0.f){
        vec3 specular_direction = 2.f*view_normal*dot(position_to_light,view_normal) - position_to_light;
        float attenuation = pow(max(dot(specular_direction,position_to_eye), 0.0),specular_falloff);
        return attenuation;
    }
    else{
        return 0.f;
    }
}

void main() {
   float diffuse_attenuation = DiffuseAttenuation();
   float specular_attenuation = SpecularAttenuation();
   vec3 composite_color = light_ambient + diffuse_attenuation*light_diffuse + specular_attenuation*light_specular;
   FragColor = vec4(composite_color,1.0);
}

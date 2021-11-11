#version 130

#define PI 3.141592

//contribution to voxel raycast:
//https://www.shadertoy.com/view/4ds3zr
//ambient occlusion and lighting from IQ
//https://www.shadertoy.com/view/4dfGzs

uniform sampler2D object_texture;

uniform vec3 params; //custom params vector (grid scale, grid extension, internal resolution)

out vec4 output_color; //new shader version

vec3 edgecolor; //check if mipmaps are working correctly

// Input from vertex shader
//in vec2 texcoord;
uniform float osg_FrameTime;



void main() {
	// pixel coordinates
	vec2 p = vec2(gl_FragCoord.x/params.z, gl_FragCoord.y / params.z);
	
	edgecolor = vec3(0.5);
	
	vec2 off = vec2(1.0 / params.z, 1.0/params.z);
	
	float vertexindex = 0.0;
		
	vec3 center = texture(object_texture, p).xyz;
	vec3 in1 = texture(object_texture, p + vec2(+off.x, 0.0)).xyz;
	
	if (in1 == center)
	{
		edgecolor.x = +off.x;
		vertexindex = vertexindex + 1.0;
	}
	
    vec3 in2 = texture(object_texture, p + vec2(-off.x, 0.0)).xyz;
	
	if (in2 == center)
	{
		edgecolor.x = -off.x;
		vertexindex = vertexindex + 1.0;
	}
	
    vec3 in3 = texture(object_texture, p + vec2(0.0, +off.y)).xyz;
	
	if (in3 == center)
	{
		edgecolor.y = +off.y;
		vertexindex = vertexindex + 1.0;
	}
	
    vec3 in4 = texture(object_texture, p + vec2(0.0, -off.y)).xyz;
	
	if (in4 == center)
	{
		edgecolor.y = -off.y;
		vertexindex = vertexindex + 1.0;
	}
	
    vec3 in5 = texture(object_texture, p + vec2(+off.x, +off.y)).xyz;
	
	if (in5 == center)
	{
		edgecolor.x = +off.x;
		edgecolor.y = +off.y;
		vertexindex = vertexindex + 1.0;
	}
	
    vec3 in6 = texture(object_texture, p + vec2(-off.x, +off.y)).xyz;
	
	if (in6 == center)
	{
		edgecolor.x = -off.x;
		edgecolor.y = +off.y;
		vertexindex = vertexindex + 1.0;
	}
	
    vec3 in7 = texture(object_texture, p + vec2(+off.x, -off.y)).xyz;
	
	if (in7 == center)
	{
		edgecolor.x = +off.x;
		edgecolor.y = -off.y;
		vertexindex = vertexindex + 1.0;
	}
	
    vec3 in8 = texture(object_texture, p + vec2(-off.x, -off.y)).xyz;
	
	if (in8 == center)
	{
		edgecolor.x = -off.x;
		edgecolor.y = -off.y;
		vertexindex = vertexindex + 1.0;
	}
	
	vec3 col = vec3(0.0, 0.0, 0.0);
	
	edgecolor = edgecolor / vertexindex;
	
	col.x = edgecolor.x * 5.0;
	col.y = edgecolor.y * 5.0;
	col.z = edgecolor.z * 5.0;
	
	output_color = vec4(col, 1.0); // final ray color + alpha channel	
}
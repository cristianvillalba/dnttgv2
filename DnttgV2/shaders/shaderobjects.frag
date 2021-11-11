#version 130

#define PI 3.141592

//contribution to voxel raycast:
//https://www.shadertoy.com/view/4ds3zr
//ambient occlusion and lighting from IQ
//https://www.shadertoy.com/view/4dfGzs

uniform sampler3D p3d_Texture0;

uniform vec3 campos; //custom campos vector
uniform vec3 target; //custom camdir vector
uniform vec3 params; //custom params vector (grid scale, grid extension, internal resolution)
uniform vec3 voxparams; //custom params vector (voxelsize, voxelsize, voxelsize)

out vec4 output_color; //new shader version

vec3 edgecolor; //check if mipmaps are working correctly

// Input from vertex shader
//in vec2 texcoord;
uniform float osg_FrameTime;

int lodvalue = 0; //3 seems a good starting point
//int lodvalue = int(osg_FrameTime) % 4;

float rand(vec2 co){
    return fract(sin(dot(co.xy ,vec2(12.9898,78.233))) * 43758.5453);
}

float sdSphere( vec3 p, float s )
{
    return length(p)-s;
}

float sdBox( vec3 p, vec3 b )
{
  vec3 q = abs(p) - b;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}



float map( in vec3 p, out vec3 hitcolor)
{
	//float noisex = rand(vec2(osg_FrameTime,0.3));
	//float noisey = rand(vec2(0.1,osg_FrameTime));
	//float noisez = rand(vec2(0.1 + osg_FrameTime,0.3));
	
	//p.x = p.x + noisex/10.0;
	//p.y = p.y + noisey/10.0;
	//p.z = p.z + noisez/10.0;
	p.y = -p.y;
	p /= (params.x*params.y);
	
	//vec4 color = texture(p3d_Texture0, p);
	vec4 color = textureLod(p3d_Texture0, p, lodvalue);
	//vec4 color = textureLod(p3d_Texture0, p, 0);
	
	hitcolor.r = color.r;
	hitcolor.g = color.g;
	hitcolor.b = color.b;
	
    return color.r;
}

// Amanatides & Woo style voxel traversal
vec3 voxelSize = vec3(params.x/voxparams.x) * (lodvalue + 1.0);//grid size divided by 10. In this case voxels of 1.0

vec2 seed;//seed for random

float simpleseed; //simple seed

vec3 worldToVoxel(vec3 i)
{
    return floor(i/voxelSize);
}

vec3 voxelToWorld(vec3 i)
{
    return i*voxelSize;	
}

vec3 voxelTrace(vec3 ro, vec3 rd, out bool hit, out vec3 hitNormal, out vec3 hitcolor)
{
    int maxSteps = 256;
	int hitindex = 0;
    vec3 voxel = worldToVoxel(ro);
    vec3 step = sign(rd);
	vec3 nearestVoxel = voxel + vec3(rd.x > 0.0, rd.y > 0.0, rd.z > 0.0);
    vec3 tMax = (voxelToWorld(nearestVoxel) - ro) / rd;
    vec3 tDelta = voxelSize / abs(rd);
    vec3 hitVoxel = voxel;
	
    hit = false;
	maxSteps = int(maxSteps / pow(2, lodvalue));
 	
    float hitT = 0.0;
    for(int i=0; i < maxSteps; i++)
	{
		if (!hit)
		{
			float d = map(voxelToWorld(voxel), hitcolor);        
			if (d != 0.0 && !hit)
			{
				hit = true;
				hitVoxel = voxel;
				hitindex = i;
                break;
			}
			bool c1 = tMax.x < tMax.y;
			bool c2 = tMax.x < tMax.z;
			bool c3 = tMax.y < tMax.z;
			if (c1 && c2) 
			{ 
				if (!hit) 
				{
					hitNormal = vec3(-step.x, 0.0, 0.0);
					hitT = tMax.x;
				}
				voxel.x += step.x;
				tMax.x += tDelta.x;
	
			} else if (c3 && !c1) 
			{
				if (!hit) 
				{
					hitNormal = vec3(0.0, -step.y, 0.0);	
					hitT = tMax.y;
				}
				voxel.y += step.y;
				tMax.y += tDelta.y;
			} else
			{
				if (!hit) 
				{
					hitNormal = vec3(0.0, 0.0, -step.z);		
					hitT = tMax.z;
				}
				voxel.z += step.z;
				tMax.z += tDelta.z;
			}
		}
    }
	
	if (hit && hitindex == maxSteps)
	{
		hit = false;
	}
	//if (hit && (hitVoxel.x > 27.0 || hitVoxel.x < -27.0 || hitVoxel.z < -27.0 || hitVoxel.z > 27.0))
	//{
	//	hit = false;
	//	return vec3(20.0);
	//}
	
	return ro + hitT * rd;
}

vec2 rand2n() {
    seed+=vec2(-1,1);
	// implementation based on: lumina.sourceforge.net/Tutorials/Noise.html
    return vec2(fract(sin(dot(seed.xy ,vec2(12.9898,78.233))) * 43758.5453),
		fract(cos(dot(seed.xy ,vec2(4.898,7.23))) * 23421.631));
};

vec2 hash2(inout float seed) {
    return fract(sin(vec2(seed+=0.1,seed+=0.1))*vec2(43758.5453123,22578.1459123));
}

vec3 ortho(vec3 v) {
    //  See : http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts
    return abs(v.x) > abs(v.z) ? vec3(-v.y, v.x, 0.0)  : vec3(0.0, -v.z, v.y);
}

vec3 getSampleBiased(vec3  dir, float power) {
	dir = normalize(dir);
	vec3 o1 = normalize(ortho(dir));
	vec3 o2 = normalize(cross(dir, o1));
	vec2 r = rand2n();
	r.x=r.x*2.*PI;
	r.y=pow(r.y,1.0/(power+1.0));
	float oneminus = sqrt(1.0-r.y*r.y);
	return cos(r.x)*oneminus*o1+sin(r.x)*oneminus*o2+r.y*dir;
}

vec3 getSample(vec3 dir) {
	return getSampleBiased(dir,0.0); // <- unbiased!
}

vec3 getCosineWeightedSample(vec3 dir) {
	return getSampleBiased(dir,1.0);
}


vec3 getBackground(vec3 dir)
{
	return vec3(0.95);
}

vec3 getConeSample(vec3 dir, float extent) {
        // Formula 34 in GI Compendium
	dir = normalize(dir);
	vec3 o1 = normalize(ortho(dir));
	vec3 o2 = normalize(cross(dir, o1));
	vec2 r =  rand2n();
	r.x=r.x*2.*PI;
	r.y=1.0-r.y*extent;
	float oneminus = sqrt(1.0-r.y*r.y);
	return cos(r.x)*oneminus*o1+sin(r.x)*oneminus*o2+r.y*dir;
}

vec3 getRayMipmap(vec3 ro, vec3 rd, out float alpha)
{
	vec3 rayOrigin = ro;
    vec3 rayDirection = rd;
	vec3 pos;
	vec3 n;
	vec3 hitcolor = vec3(0.0);
		
	lodvalue = 4; //start with a mipmap level value of 4
	
	for (int i = 4; i >= 0; i--) {
		bool hit;
		
		voxelSize = vec3((params.x/voxparams.x) * pow(2, lodvalue));
		pos = voxelTrace(rayOrigin, rayDirection, hit, n, hitcolor);
		
		if (hit){ //if hit then
			rayOrigin = pos + n *0.001* 4.0; //move rayorigin to new position
						
			if (lodvalue == 0) //if no more miplevels break
			{
				alpha = 1.0;//we will colorize this pixel
				edgecolor = hitcolor; //colorize fragment
				
						
				break;
			}
			else{
				lodvalue = lodvalue - 1; //go to the next miplevel
			}
		}
		else{
			break;
		}
	}

	return (pos + n *0.001* 4.0);
}


void main() {
	// pixel coordinates
	vec2 p = vec2((2.0*gl_FragCoord.x-params.z)/params.z, (2.0*gl_FragCoord.y-params.z)/params.z); //this will change range to [-1...1]
	
	edgecolor = vec3(0.0);
		
	float OFFSET = 0.5 * (params.x * params.y);
	
	vec3 ro = vec3( campos.x + OFFSET, campos.y -OFFSET, campos.z + OFFSET);
	vec3 ta = vec3( target.x + OFFSET, target.y -OFFSET, target.z + OFFSET);

	// camera matrix
	vec3 ww = normalize( ta - ro );
	vec3 uu = normalize( cross(ww,vec3(0.0,1.0,0.0) ) );
	vec3 vv = normalize( cross(uu,ww));
	// create view ray
	vec3 rd = normalize( p.x*uu + p.y*vv + 0.75*ww );

	//init simpleseed
	simpleseed = ro.x + ro.y * 3.43121412313;
	
	//raycast trace ray
	int RAYSAMPLES = 1;
	float alpha = 0.0;
	vec3 col = vec3(0.0);
	
	ro = getRayMipmap(ro, rd, alpha); //advance mipmap

	col = col / float(RAYSAMPLES);
	col = col +  edgecolor;//mipmap test!

	//alpha = step(0.125, col.x * col.y * col.z)*-1.0 + 1.0; //filter whites, this instead of branching 
	
	col = pow(col, vec3(0.4545)); //gamma correction
	
	output_color = vec4(col, alpha); // final ray color + alpha channel	
}
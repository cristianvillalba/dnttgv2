#version 130

#define PI 3.141592

//contribution to voxel raycast:
//https://www.shadertoy.com/view/4ds3zr
//ambient occlusion and lighting from IQ
//https://www.shadertoy.com/view/4dfGzs

uniform sampler3D p3d_Texture0;
uniform sampler2D shadow_texture;
uniform sampler2D indirect_texture;

uniform vec3 campos; //custom campos vector
uniform vec3 target; //custom camdir vector
uniform vec3 params; //custom params vector (grid scale, grid extension, internal resolution)
uniform vec3 voxparams; //custom params vector (voxelsize, voxelsize, voxelsize)

out vec4 output_color; //new shader version

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



float map( in vec3 p )
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

    return color.r;
}

// Amanatides & Woo style voxel traversal
vec3 voxelSize = vec3(params.x/voxparams.x) * (lodvalue + 1.0);//grid size divided by 10. In this case voxels of 1.0

vec2 seed;//seed for random

float simpleseed; //simple seed

vec3 mipcolor; //check if mipmaps are working correctly



vec3 worldToVoxel(vec3 i)
{
    return floor(i/voxelSize);
}

vec3 voxelToWorld(vec3 i)
{
    return i*voxelSize;	
}

vec3 voxelTrace(vec3 ro, vec3 rd, out bool hit, out vec3 hitNormal)
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
			float d = map(voxelToWorld(voxel));        
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
	
	//calc shadows here! slow has hell
	//vec3 color = vec3(1.0, 1.0, 1.0);
	//vec3 material = vec3(0.35, 0.35, 0.0);
    //float Albedo = 0.4;
	//seed = rd.xy * (osg_FrameTime + 1.0); //jitter light a little bit
	
	lodvalue = 0; //start with a mipmap level value of 4 -- doesn't work ok
	
	//cal shadows here! slow has hell
	//vec3 sun_dir = normalize(vec3(sin(osg_FrameTime*0.01), cos(osg_FrameTime*0.01), 0.0));
	
	for (int i = 0; i >= 0; i--) { // doesn't work ok with MIPMAPS
		bool hit;
		
		voxelSize = vec3((params.x/voxparams.x) * pow(2, lodvalue));
		pos = voxelTrace(rayOrigin, rayDirection, hit, n);
		
		if (hit){ //if hit then
			rayOrigin = pos + n *0.001* 4.0; //move rayorigin to new position
			//mipcolor.y += 0.05; //colorize fragment
			
			
			if (lodvalue == 0) //if no more miplevels break
			{
				alpha = 1.0;//we will colorize this pixel
				mipcolor.y = 0.5; //colorize fragment
				
				//calc shadows here! slow has hell
				//vec3 sunSampleDir = getConeSample(sun_dir,0.001);
				//float sunLight = dot(n, sunSampleDir);
				
				//voxelTrace(rayOrigin, sunSampleDir, hit, n);
											
				//if (sunLight>0.0 && !hit) {
				//	mipcolor += color*sunLight;
				//}
						
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

vec3 getRayColor(vec3 ro, vec3 rd, out float alpha, float i)
{
	vec3 color = vec3(1.0, 1.0, 1.0);
    vec3 directLight = vec3(0.0, 0.0, 0.0);
	vec3 material = vec3(0.35, 0.35, 0.0);
    float Albedo = 0.4;
    int iterations = 0;
    int object = 0;
    
    vec3 rayOrigin = ro;
    vec3 rayDirection = rd;
    
    seed = rd.xy * (osg_FrameTime + i + 1.0);
    
    
    vec3 n;
    vec3 pos;
    //vec3 sun_dir = normalize(vec3(0.0, 1.0, 0.0));
	vec3 sun_dir = normalize(vec3(sin(osg_FrameTime*0.01), cos(osg_FrameTime*0.01), 0.0));
    
    for (int i = 0; i < 1; i++) {
		bool hit;
		
		pos = voxelTrace(rayOrigin, rayDirection, hit, n);
        
        if (hit) { 
			alpha = 1.0;
			vec3 newDirection = normalize(getCosineWeightedSample(n));
            //vec3 newDirection = normalize(getSample(n));

            color *= material * Albedo * 2.0;
            
            rayOrigin = pos + n *0.001* 4.0;
            rayDirection = newDirection;
              
			vec3 sunSampleDir = getConeSample(sun_dir,0.0001);
			float sunLight = dot(n, sunSampleDir);
			
			pos = voxelTrace(rayOrigin, sunSampleDir, hit, n);
			if (sunLight>0.0 && !hit) {
				directLight += color*sunLight;
			}
            
        }
        else {
            return directLight + color * getBackground(rayDirection);
        }
    }
    
    return directLight + color * getBackground(rayDirection);
}

void main() {
	// pixel coordinates
	vec2 p = vec2((2.0*gl_FragCoord.x-params.z)/params.z, (2.0*gl_FragCoord.y-params.z)/params.z); //this will change range to [-1...1]
	vec2 pshad = vec2(gl_FragCoord.x, gl_FragCoord.y)*0.125; //reduce fragment coord to fetch pixels from other buffers
	pshad = pshad / (params.z * 0.125);//reduce fragment coord to fetch pixels from other buffers
	
	mipcolor = vec3(0.0);
		
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
	//lodvalue = 0;//reset lod params
	//voxelSize = vec3(params.x/voxparams.x);//reset voxel params
	
	//for (int i = 0; i < RAYSAMPLES; i++)
	//{	
		//vec2 rpof = 4.*(hash2(simpleseed)-vec2(0.5)) / params.z; //this will make visible gaps between billboards
		//vec2 rpof;
		//rd = normalize( (p.x+rpof.x)*uu + (p.y+rpof.y)*vv + 0.75*ww );	
		//col += getRayColor(ro, rd, alpha, float(i));
	//}

	col = col / float(RAYSAMPLES);
	col = col +  mipcolor;//mipmap test!
	//col = col + texture(shadow_texture, pshad).xyz;//apply shadow
	col = col + texture(indirect_texture, pshad).xyz*0.45;//apply indirect light
	
	//alpha = step(0.125, col.x * col.y * col.z)*-1.0 + 1.0; //filter whites, this instead of branching 
	
	col = pow(col, vec3(0.4545)); //gamma correction
	
	output_color = vec4(col, alpha); // final ray color + alpha channel	
}
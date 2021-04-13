#version 130

#define PI 3.141592
#define IQ 0
//contribution to voxel raycast:
//https://www.shadertoy.com/view/4ds3zr
//ambient occlusion and lighting from IQ
//https://www.shadertoy.com/view/4dfGzs

uniform sampler3D p3d_Texture0;

uniform vec3 campos; //custom campos vector
uniform vec3 target; //custom camdir vector
uniform vec3 params; //custom params vector (scale, scale, internal resolution)

// Input from vertex shader
//in vec2 texcoord;
uniform float osg_FrameTime;

float rand(vec2 co){
    return fract(sin(dot(co.xy ,vec2(12.9898,78.233))) * 43758.5453);
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
	p *= (1/params.x);
	

	vec4 color = texture(p3d_Texture0, p);
		
    return color.r;
}

float sdSphere( vec3 p, float s )
{
    return length(p)-s;
}

// Amanatides & Woo style voxel traversal
//vec3 voxelSize = vec3(sin(osg_FrameTime*0.25) + 1.0); // in world space
vec3 voxelSize = vec3(0.05*params.x);//original value

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

vec3 voxelTrace(vec3 ro, vec3 rd, out bool hit, out vec3 hitNormal, out vec3 outvox)
{
    const int maxSteps = 64;
    const float isoValue = 0.0;
	
	float rvalue = 0.0;
	float gvalue = 0.0;
	float bvalue = 0.0;

    vec3 voxel = worldToVoxel(ro);
    vec3 step = sign(rd);

    vec3 nearestVoxel = voxel + vec3(rd.x > 0.0, rd.y > 0.0, rd.z > 0.0);
    vec3 tMax = (voxelToWorld(nearestVoxel) - ro) / rd;
    vec3 tDelta = voxelSize / abs(rd);

    vec3 hitVoxel = voxel;
	
    hit = false;
    float hitT = 0.0;
    for(int i=0; i<maxSteps; i++) {
        float d = map(voxelToWorld(voxel));        
        if (d != isoValue && !hit) {
            hit = true;
	    	hitVoxel = voxel;
			outvox = hitVoxel;
            break;
        }

        if (tMax.x < tMax.y && tMax.x < tMax.z) { 
            voxel.x += step.x;
            tMax.x += tDelta.x;
			if (!hit) {
				
				if (-step.x == -1){
					rvalue = -1.0;
					gvalue = 0.0;
					bvalue = 0.0;
				}
				else 
				{
					rvalue = 1.0;
					gvalue = 0.0;
					bvalue = 0.0;
				}
				
				hitT = tMax.x;
			}
        } else if (tMax.y < tMax.z) {
            voxel.y += step.y;
            tMax.y += tDelta.y;
			if (!hit) {
				
				if (-step.y == -1){
					rvalue = 0.0;
					gvalue = -1.0;
					bvalue = 0.0;
				}
				else 
				{
					rvalue = 0.0;
					gvalue = 1.0;
					bvalue = 0.0;
				}
				
				hitT = tMax.y;
			}
        } else {
            voxel.z += step.z;
            tMax.z += tDelta.z;
			if (!hit) {
				
				if (-step.z == -1){
					rvalue = 0.0;
					gvalue = 0.0;
					bvalue = -1.0;
				}
				else 
				{
					rvalue = 0.0;
					gvalue = 0.0;
					bvalue = 1.0;
				}
				
				hitT = tMax.z;
			}
        }
         
    }
	
	hitNormal = vec3(rvalue, gvalue, bvalue); //with black colors

	//outvox = voxelToWorld(hitVoxel);
	outvox = hitVoxel;
    //return voxelToWorld(hitVoxel);
	return ro + hitT*rd;
}

bool lightTrace(vec3 ro, vec3 rd)
{
	float t = 0;
	
    for (int i = 0; i < 200; i++)
	{
		vec3 pos = ro + t * rd;
		
		float distance = map(pos);
		
		if (distance != 0.0 && i > 20) //prevent hitting the same object(?)
		{
			return true;
		}
		
		t += 0.1;
	}
	
	return false;
}

float calcOcc( in vec2 uv, vec4 va, vec4 vb, vec4 vc, vec4 vd )
{
    vec2 st = 1.0 - uv;

    // edges
    vec4 wa = vec4( uv.x, st.x, uv.y, st.y ) * vc;

    // corners
    vec4 wb = vec4(uv.x*uv.y,
                   st.x*uv.y,
                   st.x*st.y,
                   uv.x*st.y)*vd*(1.0-vc.xzyw)*(1.0-vc.zywx);
    
    return wa.x + wa.y + wa.z + wa.w +
           wb.x + wb.y + wb.z + wb.w;
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
	return vec3(1.0);
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

vec3 getRayColor(vec3 ro, vec3 rd, out float alpha)
{
	vec3 luminance = vec3(1.5);
	int RayDepth = 2;
    bool hit;

    vec3 n; //normal
	vec3 outvox; //hit out voxel
	
	seed = rd.xy + vec2(osg_FrameTime + 1.0,osg_FrameTime - 1.0);
	
	vec3 material = vec3(3.0, 2.0, 1.5);
	vec3 direct = vec3(0.0);
	//vec3 sun_dir = normalize(vec3(sin(osg_FrameTime*0.5), cos(osg_FrameTime*0.5),0.0));
	vec3 sun_dir = normalize(vec3(1.0, 1.0, 0.0));
	
	
	//Just render normals
	vec3 pos = voxelTrace(ro, rd, hit, n, outvox);
		
	if (hit) {
		alpha = 1.0;
		
		return (n+vec3(1.,1.,1.))/2.; //with all colors;
	}

	return vec3(0.0);
	
	 
	/* for (int i=0; i < RayDepth; i++) {
		vec3 pos = voxelTrace(ro, rd, hit, n, outvox);
		
		if (hit) {
			alpha = 1.0;
			//rd = getSample(n); // new direction (towards light)
			rd = getCosineWeightedSample(n); // new direction (towards light)
			luminance *= material*0.125; //0.125 is albedo
			ro = voxelToWorld(outvox) + n*1.0001; // new start point
			
			// Direct lighting
			vec3 sunSampleDir = getConeSample(sun_dir,1E-5);
			float sunLight = dot(n, sunSampleDir);
			
			pos = voxelTrace(voxelToWorld(outvox) + n*1.0001, sunSampleDir, hit, n, outvox);
			if (sunLight>0.0 && !hit) {
				direct += luminance*sunLight;
			}
		}else{
			return direct + luminance * getBackground( rd );
		}
	}
	
	return direct + luminance * getBackground( rd ); */
}

void main() {
  
	// pixel coordinates
	vec2 p = vec2((2.0*gl_FragCoord.x-params.z)/params.z, (2.0*gl_FragCoord.y-params.z)/params.z); //this will change range to [-1...1]
	
	// camera movement	
	float an = 0.2*osg_FrameTime;
	
	float OFFSET = 0.5 * params.x;
	
	vec3 ro = vec3( campos.x + OFFSET, campos.y -OFFSET, campos.z + OFFSET);
	vec3 ta = vec3( target.x + OFFSET, target.y -OFFSET, target.z + OFFSET);

	// camera matrix
	vec3 ww = normalize( ta - ro );
	vec3 uu = normalize( cross(ww,vec3(0.0,1.0,0.0) ) );
	vec3 vv = normalize( cross(uu,ww));
	// create view ray
	vec3 rd = normalize( p.x*uu + p.y*vv + 0.75*ww );

	if (IQ == 0)
	{
		//init simpleseed
		simpleseed = ro.x + ro.y * 3.43121412313;
		
		//raycast trace ray
		int RAYSAMPLES = 1;
		float alpha = 0.0;
		vec3 col = vec3(0.0);
		
		for (int i = 0; i < RAYSAMPLES; i++)
		{	
			//vec2 rpof = 4.*(hash2(simpleseed)-vec2(0.5)) / params.z; //this will make some edges visible in billboards
			vec2 rpof;
			rd = normalize( (p.x+rpof.x)*uu + (p.y+rpof.y)*vv + 0.75*ww );
			col += getRayColor(ro, rd, alpha);
		}

		col = col / float(RAYSAMPLES);	
		
		gl_FragColor = vec4(col, alpha); // final ray color + alpha channel
	}
	else
	{
		//IQ lighting 
		bool hit;

		vec3 n; //normal
		vec3 outvox; //hit out voxel
		vec3 col = vec3(0.0);
	
		vec3 pos = voxelTrace(ro, rd, hit, n, outvox);
		
		if (hit)
		{		
			vec3 uvw = pos - voxelToWorld(outvox);
			
			vec3 material = vec3(0.2);
			
			//vec3 sun_dir = normalize(vec3(1.0, 1.0, -1.0));
			vec3 newn = n * vec3(-1.0, 1.0, -1.0); //flip y;
			//vec3 sun_dir = normalize(vec3(0.0,sin(osg_FrameTime*0.5), cos(osg_FrameTime*0.5)));
			vec3 sun_dir = normalize(vec3(sin(osg_FrameTime*0.5), cos(osg_FrameTime*0.5), sin(osg_FrameTime*0.5)));
			float sun_diff = clamp( dot(newn, sun_dir), 0.0, 1.0); //dot product with sun and normal
			float sky_diff = clamp( 0.5 + 0.5*dot(newn, vec3(0.0, 1.0, 0.0)), 0.0, 1.0); //dot product with sky(like a light coming from Y axis) and normal + bias --- change from -1 - 1 to 0 - 1 
			float bounce_diff = clamp( 0.5 + 0.5*dot(newn, vec3(0.0, -1.0, 0.0)), 0.0, 1.0);

			vec3 newsundir = sun_dir * vec3(-1.0, 1.0, -1.0);
			
			bool hitlight = lightTrace(pos + newn*0.1, newsundir);
			float sun_shad = 1.0;
			
			if (hitlight)
			{
				sun_shad = 0.0;
			}
				
			vec3 v1  = voxelToWorld(outvox + n + n.yzx);
			vec3 v2  = voxelToWorld(outvox + n - n.yzx);
			vec3 v3  = voxelToWorld(outvox + n + n.zxy);
			vec3 v4  = voxelToWorld(outvox + n - n.zxy);
			vec3 v5  = voxelToWorld(outvox + n + n.yzx + n.zxy);
			vec3 v6  = voxelToWorld(outvox + n - n.yzx + n.zxy);
			vec3 v7  = voxelToWorld(outvox + n - n.yzx - n.zxy);
			vec3 v8  = voxelToWorld(outvox + n + n.yzx - n.zxy);
			vec3 v9  = voxelToWorld(outvox + n.yzx);
			vec3 v10 = voxelToWorld(outvox - n.yzx);
			vec3 v11 = voxelToWorld(outvox + n.zxy);
			vec3 v12 = voxelToWorld(outvox - n.zxy);
			vec3 v13 = voxelToWorld(outvox + n.yzx + n.zxy); 
			vec3 v14 = voxelToWorld(outvox - n.yzx + n.zxy);
			vec3 v15 = voxelToWorld(outvox - n.yzx - n.zxy);
			vec3 v16 = voxelToWorld(outvox + n.yzx - n.zxy);

			vec4 vc = vec4( map(v1),  map(v2),  map(v3),  map(v4)  );
			vec4 vd = vec4( map(v5),  map(v6),  map(v7),  map(v8)  );
			vec4 va = vec4( map(v9),  map(v10), map(v11), map(v12) );
			vec4 vb = vec4( map(v13), map(v14), map(v15), map(v16) );
			
			vec2 uv = vec2( dot(n.yzx, uvw), dot(n.zxy, uvw) );
			
			float occ = 1.0;
			
			// ambient occlusion
			occ = calcOcc( uv, va, vb, vc, vd );
			occ = 1.0 - occ/16.0;
			occ = occ*occ;
			occ = occ*occ;
			
			col = material*vec3(3.0, 2.0, 1.5) * sun_diff*sun_shad;
			col += material*vec3(0.45, 0.86, 0.89) * sky_diff;
			col += material*vec3(0.7, 0.3, 0.2) * bounce_diff;
			
			col *= occ;//ambient occlusion
			
			col = pow(col, vec3(0.4545)); //gamma correction
			
			//col = vec3(occ);
			//col = pos/params.x;
			
			gl_FragColor = vec4(col, 1.0);
		}
		else
		{
			gl_FragColor = vec4(1, 1, 1, 0);
		} 
		
		//* IQ lighting 
	}	
}
#version 130

#define PI 3.141592
//https://www.shadertoy.com/view/4ds3zr

//#define SCALE 4.0
//#define OFFSET 0.5 * SCALE //it was 0.5

uniform sampler3D p3d_Texture0;

uniform vec3 campos; //custom campos vector
uniform vec3 target; //custom camdir vector
uniform vec3 scale;  //custom scale  vector

// Input from vertex shader
//in vec2 texcoord;
uniform float osg_FrameTime;

mat4 scalematrix = mat4(5.0, 0.0, 0.0, 0.0,  // 1. column
						0.0, 5.0, 0.0, 0.0,  // 2. column
						0.0, 0.0, 5.0, 0.0,  // 3. column
						0.0, 0.0, 0.0, 1.0); // 4. column

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
	p *= (1/scale.x);
	

	vec4 color = texture(p3d_Texture0, p);
	
    return color.r;
}

float sdSphere( vec3 p, float s )
{
    return length(p)-s;
}


// lighting
vec3 shade(vec3 pos, vec3 n, vec3 eyePos)
{
    const vec3 lightPos = vec3(4.0, 3.0, 5.0);
    const vec3 color = vec3(1.0, 1.0, 0.0);
    const float shininess = 40.0;

    vec3 l = normalize(lightPos - pos);
    vec3 v = normalize(eyePos - pos);
    vec3 h = normalize(v + l);
    float diff = dot(n, l);
    float spec = max(0.0, pow(dot(n, h), shininess)) * float(diff > 0.0);
    diff = max(0.0, diff);
    //diff = 0.5+0.5*diff;

    float fresnel = pow(1.0 - dot(n, v), 5.0);
    //float ao = ambientOcclusion(pos, n);

	return vec3(diff);
//    return vec3(diff*ao)*color + vec3(spec + fresnel*0.5);
//    return vec3(diff*ao)*color;	
//    return vec3(diff*ao)*color + vec3(spec);
//    return vec3(ao);
//    return vec3(fresnel);
}

// Amanatides & Woo style voxel traversal
//vec3 voxelSize = vec3(sin(osg_FrameTime) + 1.0); // in world space
vec3 voxelSize = vec3(0.05*scale.x);


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
    const int maxSteps = 64;
    const float isoValue = 0.0;
	
	float rvalue = 0;
	float gvalue = 0;
	float bvalue = 0.5;

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
            //break;
        }

        if (tMax.x < tMax.y && tMax.x < tMax.z) { 
            voxel.x += step.x;
            tMax.x += tDelta.x;
			if (!hit) {
				//hitNormal = vec3(-step.x, 0.0, 0.0);
				
				if (-step.x == -1){
					//hitNormal = vec3(1.0, 0.0, 0.0);
					rvalue = 1.0;
					gvalue = 0.5;
					bvalue = 0.5;
				}
				else 
				{
					//hitNormal = vec3(0.0, 0.0, 0.0);
					rvalue = 0.0;
					gvalue = 0.5;
					bvalue = 0.5;
				}
				
				hitT = tMax.x;
			}
        } else if (tMax.y < tMax.z) {
            voxel.y += step.y;
            tMax.y += tDelta.y;
			if (!hit) {
				//hitNormal = vec3(0.0, -step.y, 0.0);		
				//hitNormal = vec3(0.0, 1.0, 0.0);		
				
				if (-step.y == -1){
					//hitNormal = vec3(1.0, 0.0, 0.0);
					rvalue = 0.5;
					gvalue = 1.0;
					bvalue = 0.5;
				}
				else 
				{
					//hitNormal = vec3(0.0, 0.0, 0.0);
					rvalue = 0.5;
					gvalue = 0.0;
					bvalue = 0.5;
				}
				
				hitT = tMax.y;
			}
        } else {
            voxel.z += step.z;
            tMax.z += tDelta.z;
			if (!hit) {
				//hitNormal = vec3(0.0, 0.0, -step.z);
				//hitNormal = vec3(0.0, 0.0, 1.0);
				
				if (-step.z == -1){
					//hitNormal = vec3(1.0, 0.0, 0.0);
					rvalue = 0.5;
					gvalue = 0.5;
					bvalue = 1.0;
				}
				else 
				{
					//hitNormal = vec3(0.0, 0.0, 0.0);
					rvalue = 0.5;
					gvalue = 0.5;
					bvalue = 1.0;
				}
				
				hitT = tMax.z;
			}
        }
         
    }
	
	hitNormal = vec3(rvalue, gvalue, bvalue);

    //return voxelToWorld(hitVoxel);
	return ro + hitT*rd;
}

void main() {
	//vec4 color = texture(p3d_Texture0, texcoord);
	  
	// pixel coordinates
	//vec2 p = vec2((2.0*gl_FragCoord.x-800)/600, (2.0*gl_FragCoord.y-600)/600);
	vec2 p = vec2((2.0*gl_FragCoord.x-128)/128, (2.0*gl_FragCoord.y-128)/128); //this will change range to [-1...1]
	
	// camera movement	
	float an = 0.2*osg_FrameTime;
	
	float OFFSET = 0.5 * scale.x;
	//float OFFSET = 0.0;
	
	//center of texture3d is 0.5, 0.5, 0.5
	//vec3 ro = vec3( 1.0*sin(an) + 0.5, -0.5, 1.0*cos(an) + 0.5); //this will spin in center
	//vec3 ro = vec3( 0.5, -0.5, 1.5); //target of the camera	
	vec3 ro = vec3( campos.x + OFFSET, campos.y -OFFSET, campos.z + OFFSET);
	//vec3 ta = vec3( 0.5, -0.5, 0.5 );
	vec3 ta = vec3( target.x + OFFSET, target.y -OFFSET, target.z + OFFSET);

	// camera matrix
	vec3 ww = normalize( ta - ro );
	vec3 uu = normalize( cross(ww,vec3(0.0,1.0,0.0) ) );
	vec3 vv = normalize( cross(uu,ww));
	// create view ray
	vec3 rd = normalize( p.x*uu + p.y*vv + 1.5*ww );

	//raycast
	// trace ray
    bool hit;

    vec3 n;
    vec3 pos = voxelTrace(ro, rd, hit, n);
	
	if (hit)
	{
		// shade
        vec3 rgb = shade(pos, n, ro);
		//gl_FragColor = vec4(0.8, 0, 0, 1.0);
		//gl_FragColor = vec4(rgb.r, rgb.g, rgb.b, 1.0); //some lighting
		gl_FragColor = vec4(n.x, n.y, n.z, 1.0); //check normals
	}
	else
	{
		gl_FragColor = vec4(1, 1, 1, 0);
	}

	
}
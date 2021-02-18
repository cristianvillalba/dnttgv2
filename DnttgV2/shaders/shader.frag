#version 130

#define PI 3.141592

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
	
	//vec4 prevpos = vec4( p.x, p.y , p.z, 1 );
	//vec4 newpos = prevpos * scalematrix;
	//vec3 finalpos = vec3(newpos.x, newpos.y, newpos.z);
		
	//vec4 color = texture(p3d_Texture0, finalpos);
	vec4 color = texture(p3d_Texture0, p);
	
    return color.r;
}

float sdSphere( vec3 p, float s )
{
    return length(p)-s;
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

	// raymarch
	const float tmax = 128;
	float t = 0.0;
	for( int i=0; i<128; i++ )
	{
		vec3 pos = ro + rd*40*i/128;
		pos.y = -pos.y;
		pos *= 1/scale.x;
		
		//pos *= 1/(sin(osg_FrameTime)*1.0 + 1.0);
	
		//pos *= 1/scale.x;
		float h = map(pos);
		//float h = sdSphere(pos, 0.25);
		if( h != 0.0 ){ 
			//t = (1.0 - i/100.0);
			//break;
			//t = (1.0 - i/150.0);
			t = (1.0 - i/75.0);
			//t = h;
			break;
		}
	}
  
	//vec2 pos = vec2((gl_FragCoord.x/800*2.0 -1.0) * PI, (gl_FragCoord.y/600*2.0 - 1.0)* PI);
	//gl_FragColor = vec4(0.2, 0.6, 1., 1.) * abs(sin(20.*pos.y + 20.*sin(pos.x + osg_FrameTime)));
		
	//if t != 0 then alpha is 1.0
	//otherwise alpha is 0.0
	gl_FragColor = vec4(t, 0, 0, abs(step(0,-t) - 1));
	
}
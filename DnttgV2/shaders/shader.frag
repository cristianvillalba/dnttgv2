#version 130

#define PI 3.141592

uniform sampler3D p3d_Texture0;

uniform vec3 campos; //custom campos vector
uniform vec3 target; //custom camdir vector

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
	vec2 p = vec2((2.0*gl_FragCoord.x-800)/600, (2.0*gl_FragCoord.y-600)/600);
	
	// camera movement	
	float an = 0.2*osg_FrameTime;
	
	//center of texture3d is 0.5, 0.5, 0.5
	//vec3 ro = vec3( 1.0*sin(an) + 0.5, -0.5, 1.0*cos(an) + 0.5); //this will spin in center
	//vec3 ro = vec3( 0.5, -0.5, 1.5); //target of the camera	
	vec3 ro = vec3( campos.x + 0.5, campos.y -0.5, campos.z + 0.5);
	//vec3 ta = vec3( 0.5, -0.5, 0.5 );
	vec3 ta = vec3( target.x + 0.5, target.y -0.5, target.z + 0.5);
	
	// camera matrix
	vec3 ww = normalize( ta - ro );
	vec3 uu = normalize( cross(ww,vec3(0.0,1.0,0.0) ) );
	vec3 vv = normalize( cross(uu,ww));
	// create view ray
	vec3 rd = normalize( p.x*uu + p.y*vv + 1.5*ww );

	// raymarch
	const float tmax = 128;
	float t = 0.0;
	for( int i=0; i<64; i++ )
	{
		vec3 pos = ro + rd*2*i/64;
		pos.y = -pos.y;
		float h = map(pos);
		//float h = sdSphere(pos, 0.25);
		if( h != 0.0 ){ 
			//t = (1.0 - i/100.0);
			//break;
			t = (1.0 - i/150.0);
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
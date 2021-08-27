// DnttgV2.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <map>

#include "pandaFramework.h"
#include "pandaSystem.h"
#include "geomPoints.h"
#include "geomTriangles.h"

#include "shader.h"

#include "genericAsyncTask.h"
#include "asyncTaskManager.h"

#include "cIntervalManager.h"
#include "cLerpNodePathInterval.h"
#include "cMetaInterval.h"

#include <load_prc_file.h>

#include "orthographicLens.h"
#include <texturePool.h>

#include "VDBGrid.h"
#include "FilterManager.h"

#include "antialiasAttrib.h"

#include "glgsg.h"
#include "gl/GL.h"
#include "gl/GLU.h"
#include "gl/glext.h"



#include "Dnntgv2.h"


#define WIDTH 200  //for bunny - for first raycaster
#define HEIGHT 200 //for bunny - for first raycaster
#define INTERNALRES 256 //internal texture resolution
#define BUNNY 1 //old vs new raycaster
#define BOUNDINGBOX 1 //bounding box of 3d texture
#define TEXTURESIZE 32 //3d texture resolution
#define GRIDEXTENSION 1 //how many side voxels this will render
#define DENOISE 1 //denoising shader as an image post processing

using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::milli;
using std::random_device;
using std::sort;
using std::vector;

// The global task manager
PT(AsyncTaskManager) taskMgr = AsyncTaskManager::get_global_ptr();
// The global clock
PT(ClockObject) globalClock = ClockObject::get_global_clock();
// Here's what we'll store the camera in.
NodePath camera;

//Screen Projection Plane - for rasterizer with points
PT(GeomVertexData) vdata;

//VDB class handler
VDBGrid*  grid;

//global 3d texture
PT(Texture) bunn;

//Empty 3d texture;
PT(Texture) emptyTexture;
unsigned char * emptyTextureArray;

NodePath mainQuad[1];
NodePath mainQuadNorm[1];

unsigned char * gridTextureArray[((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)];
unsigned char * gridTextureArrayBuffer[((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)];
UsedTextures usedTexturesVector;

LVector3f gridCoordOffset[((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)];

float offsetvectorx[((GRIDEXTENSION * 2) + 1)];
float offsetvectory[((GRIDEXTENSION * 2) + 1)];
float offsetvectorz[((GRIDEXTENSION * 2) + 1)];

float offsetvectorxO[((GRIDEXTENSION * 2) + 1)];
float offsetvectoryO[((GRIDEXTENSION * 2) + 1)];
float offsetvectorzO[((GRIDEXTENSION * 2) + 1)];

int gridoffsetx[((GRIDEXTENSION * 2) + 1)];
int gridoffsety[((GRIDEXTENSION * 2) + 1)];
int gridoffsetz[((GRIDEXTENSION * 2) + 1)];

//Window main
WindowFramework *mainWindow;

//Framework main
PandaFramework *mainFramework;

//Global camera pos;
float CAM_x = 0.0;
float CAM_y = 0.0;
float CAM_z = 0.0;

//Global camera pos;
int GRID_x = 0;
int GRID_y = 0;
int GRID_z = 0;

//Grid parameters
float GRID_SCALE = 10.0f;
float TEXTURE_3D_EXTENSION = (GRIDEXTENSION * 2.0) + 1.0;
float VOXEL_SIZE = 20.0f; //this number will divide the GRID_SCALE to get the actual voxel size

//grid frustrum
GridFrustrum gridFrustrum;

//Velocity vector/Acceleration
float VELX = 0.0;
float VELY = 0.0;
float VELZ = 0.0;
float ACCELERATION = 20.0;
float FRICTION = 0.2;
float SPINVEL = 0;
float ANGLEDEGREES = 0;

//key states
bool FORWARD = false;
bool BACKWARDS = false;
bool LEFT = false;
bool RIGHT = false;
bool UP = false;
bool DOWN = false;
bool SPINL = false;
bool SPINR = false;
bool SPAWN = false;

//refresh flag
bool REFRESHGRID = false;

//Mutex to access critical section
Mutex mutex("CriticalStuff");
bool* pREFRESHGRID = &REFRESHGRID;

//Chain to parallel refreshing
AsyncTaskChain * renderChain;

std::queue<KeyTriple> renderQueue;
std::queue<KeyTriple> refreshQueue;


//---OpenGLExtensions
void *GetAnyGLFuncAddress(const char *name)
{
	void *p = (void *)wglGetProcAddress(name);
	if (p == 0 ||
		(p == (void*)0x1) || (p == (void*)0x2) || (p == (void*)0x3) ||
		(p == (void*)-1))
	{
		HMODULE module = LoadLibraryA("opengl32.dll");
		p = (void *)GetProcAddress(module, name);
	}

	return p;
}

PFNGLTEXTURESUBIMAGE3DPROC glTextureSubImage3D;
//PFNGLTEXTURESUBIMAGE3DPROC glTextureSubImage3D = (PFNGLTEXTURESUBIMAGE3DPROC)GetAnyGLFuncAddress("glTextureSubImage3D");


PT(Texture) rendertexture01;
PT(Texture) rendertexture02;
PT(Texture) rendertexture03;
PT(Texture) prevframe_texture;

NodePath tempAAbillboard; //temporal AA shader billboard


void print_results(const char *const tag,
	high_resolution_clock::time_point startTime,
	high_resolution_clock::time_point endTime) {
	printf("%s: Time: %fms\n", tag, duration_cast<duration<double, milli>>(endTime - startTime).count());
}


void translate(float tx, float ty, float tz, float * txout, float * tyout, float * tzout) {
	*txout += tx;
	*tyout += ty;
	*tzout += tz;
}


void rotatex(float angle, float * txout, float * tyout, float * tzout) {
	angle = angle * M_PI / 180.0;
	float temp = *tyout;
	*tyout = *tyout * cos(angle) - *tzout * sin(angle);
	*tzout = temp * sin(angle) + *tzout * cos(angle);
}

void rotatey(float angle, float * txout, float * tyout, float * tzout) {
	angle = (angle * M_PI) / 180.0;
	float temp = *tzout;
	*tzout = *tzout * cos(angle) - *txout * sin(angle);
	*txout = temp * sin(angle) + *txout * cos(angle);
}

void rotatez(float angle, float * txout, float * tyout, float * tzout) {
	angle = angle * M_PI / 180.0;
	float temp = *txout;
	*txout = *txout * cos(angle) - *tyout * sin(angle);
	*tyout = temp * sin(angle) + *tyout *cos(angle);
}

void scale(float sf, float xf, float yf, float zf, float * txout, float * tyout, float * tzout) {
	*txout = *txout * sf + (1 - sf) * xf;
	*tyout = *tyout * sf + (1 - sf) * yf;
	*tzout = *tzout * sf + (1 - sf) * zf;
}


void advanceCamera(const Event* eventPtr, void* dataPtr)
{
	if (!FORWARD)
	{
		FORWARD = true;
	}
}

void advanceCameraUp(const Event* eventPtr, void* dataPtr)
{
	if (FORWARD)
	{
		FORWARD = false;
	}
}

void backCamera(const Event* eventPtr, void* dataPtr)
{
	if (!BACKWARDS)
	{
		BACKWARDS = true;
	}

}

void backCameraUp(const Event* eventPtr, void* dataPtr)
{
	if (BACKWARDS)
	{
		BACKWARDS = false;
	}
}

void leftCamera(const Event* eventPtr, void* dataPtr)
{
	if (!LEFT)
	{
		LEFT = true;
	}
}

void leftCameraUp(const Event* eventPtr, void* dataPtr)
{
	if (LEFT)
	{
		LEFT = false;
	}
}

void rightCamera(const Event* eventPtr, void* dataPtr)
{
	if (!RIGHT)
	{
		RIGHT = true;
	}
}

void rightCameraUp(const Event* eventPtr, void* dataPtr)
{
	if (RIGHT)
	{
		RIGHT = false;
	}
}

void upCamera(const Event* eventPtr, void* dataPtr)
{
	if (!UP)
	{
		UP = true;
	}
}

void upCameraUp(const Event* eventPtr, void* dataPtr)
{
	if (UP)
	{
		UP = false;
	}
}

void downCamera(const Event* eventPtr, void* dataPtr)
{
	if (!DOWN)
	{
		DOWN = true;
	}
}

void downCameraUp(const Event* eventPtr, void* dataPtr)
{
	if (DOWN)
	{
		DOWN = false;
	}
}

void spinLCamera(const Event* eventPtr, void* dataPtr)
{
	if (!SPINL)
	{
		SPINL = true;
	}
}

void spinLCameraUp(const Event* eventPtr, void* dataPtr)
{
	if (SPINL)
	{
		SPINL = false;
	}
}

void spinRCamera(const Event* eventPtr, void* dataPtr)
{
	if (!SPINR)
	{
		SPINR = true;
	}
}

void spinRCameraUp(const Event* eventPtr, void* dataPtr)
{
	if (SPINR)
	{
		SPINR = false;
	}
}

void spawnSphere(const Event* eventPtr, void* dataPtr)
{
	if (!SPAWN)
	{
		SPAWN = true;
	}
}


AsyncTask::DoneStatus modifyGrid(GenericAsyncTask *task, void *data)
{
	if (SPAWN)
	{
		LVector3f fw = mainWindow->get_render().get_relative_point(camera, LVector3f(0, 0, 5.0));
		//LVector3f fw = mainWindow->get_render().get_relative_point(camera, LVector3f(0, 0, 0));
		fw /= GRID_SCALE;
		fw *= -1000;
		fw.add_x(GRID_x * 1000);
		fw.add_y(GRID_y * 1000);
		fw.add_z(GRID_z * 1000);

		std::cout << "spawn at x: " << fw.get_x() << "spawn at y: " << fw.get_y() << "spawn at z: " << fw.get_z() << "\n";

		fw *= -1; //invert coordinates to spawn
		grid->spawnSphere(fw, 200.0f);

		//only refresh center
		KeyTriple params = std::make_tuple(GRID_x, GRID_y, GRID_z);
		refresh3dTextureAsArray(gridFrustrum[params], GRID_x, GRID_y, GRID_z);
		refreshQueue.push(params);

		SPAWN = false;
	}

	return AsyncTask::DS_cont;
}

AsyncTask::DoneStatus refreshGrid(GenericAsyncTask *task, void *data)
{
	while (refreshQueue.size() > 0) {
		KeyTriple refresh = refreshQueue.front();
		refreshQueue.pop();

		RefreshTexture(refresh);
	}

	return AsyncTask::DS_cont;
}

PT(Texture) Render3dTexture(int gridx, int gridy, int gridz)
{
	int texsize = TEXTURESIZE;

	PT(Texture)  bunn = new Texture("bunn");
	bunn->setup_3d_texture(texsize, texsize, texsize, Texture::ComponentType::T_float, Texture::Format::F_rgba8);

	//std::cout << " grid x: " << gridx << " y: " << gridy << " z: " << gridz << "  \n";

	for (int k = 0; k < texsize; k++) {
		PNMImage* pPNMImage = new PNMImage(texsize, texsize, 4);

		for (int i = 0; i < texsize; i++)
		{
			for (int j = 0; j < texsize; j++)
			{
				//This is good for index sample
				float x = (((float)i / texsize) - 0.5f) * 1000.0f;
				float y = (((float)j / texsize) - 0.5f) * 1000.0f; //250 to put the bunny in the center of the render
				float z = (((float)k / texsize) - 0.5f) * 1000.0f;

				//this is good for world coordinates sample
				//float x = (((float)i / texsize) - 0.5f) * 100.0f;
				//float y = (((float)j / texsize) - 0.5f) * 100.0f + 25; //250 to put the bunny in the center of the render
				//float z = (((float)k / texsize) - 0.5f) * 100.0f;

				float r = 0.0f;
				float g = 0.0f;
				float b = 0.0f;

				//translate(500.0, 0, 0, &x, &y, &z);
				//float data = grid->getValue(x, y, z);
				float data = grid->getValue(x - gridx * 1000, y - gridy * 1000, z - gridz * 1000);

				if (data > 0)
				{
					r = 1.0f;
					g = 1.0f;
				}

				if (BOUNDINGBOX == 1)
				{
					if ((i == 0 || i == (texsize - 1)) && (k == 0 || k == (texsize - 1)) ||
						(j == 0 || j == (texsize - 1)) && (k == 0 || k == (texsize - 1)) ||
						(i == 0 && j == 0) || (i == 0 && j == (texsize - 1)) ||
						(i == (texsize - 1) && j == 0) || (i == (texsize - 1) && j == (texsize - 1))
						)
					{
						r = 1.0f;
						g = 1.0f;
					}
				}

				//int flipy = -j + texsize - 1;
				//std::cout << flipy << "\n";

				pPNMImage->set_red(i, j, r);
				pPNMImage->set_green(i, j, g);
				pPNMImage->set_blue(i, j, b);
				pPNMImage->set_alpha(i, j, 1.0f);
			}
		}

		bunn->load(*pPNMImage, k, 0);

		delete pPNMImage;
	}

	//bunn->write(Filename("woodgrain-#.png"), 0, 0, true, false);

	bunn->set_wrap_u(SamplerState::WrapMode::WM_border_color);
	bunn->set_wrap_v(SamplerState::WrapMode::WM_border_color);
	bunn->set_wrap_w(SamplerState::WrapMode::WM_border_color);
	bunn->set_border_color(LColor(0.0, 0.0, 0.0, 0.0));
	bunn->set_keep_ram_image(true);
	return bunn;
}

unsigned char * Render3dTextureAsArray(int gridx, int gridy, int gridz)
{
	int texsize = TEXTURESIZE;

	unsigned char * tex;
	tex = (unsigned char *)malloc(texsize * texsize * texsize * 4 * (sizeof(char)));

	//std::cout << " grid x: " << gridx << " y: " << gridy << " z: " << gridz << "  \n";
	int n = 0;

	//for (int k = 0; k < texsize; k++) {

	//	for (int i = 0; i < texsize; i++)
	//	{
	//		for (int j = 0; j < texsize; j++)
	//		{
	//			//This is good for index sample
	//			float x = (((float)i / texsize) - 0.5f) * -1000.0f; //invert axis
	//			float y = (((float)j / texsize) - 0.5f) * 1000.0f;
	//			float z = (((float)k / texsize) - 0.5f) * 1000.0f;

	//			//this is good for world coordinates sample
	//			//float x = (((float)i / texsize) - 0.5f) * 100.0f;
	//			//float y = (((float)j / texsize) - 0.5f) * 100.0f + 25; //250 to put the bunny in the center of the render
	//			//float z = (((float)k / texsize) - 0.5f) * 100.0f;

	//			int r = 0;
	//			int g = 0;
	//			int b = 0;

	//			float data = grid->getValue(y - gridy * 1000, x - gridx * 1000, z - gridz * 1000);

	//			//if (data > 0)
	//			//std::cout << "value: " << data << "\n";
	//			if (data != 1.5) //why 1.5?
	//			{
	//				r = 255;
	//				g = 255;
	//			}

	//			if (BOUNDINGBOX == 1)
	//			{
	//				if ((i == 0 || i == (texsize - 1)) && (k == 0 || k == (texsize - 1)) ||
	//					(j == 0 || j == (texsize - 1)) && (k == 0 || k == (texsize - 1)) ||
	//					(i == 0 && j == 0) || (i == 0 && j == (texsize - 1)) ||
	//					(i == (texsize - 1) && j == 0) || (i == (texsize - 1) && j == (texsize - 1))
	//					)
	//				{
	//					r = 255;
	//					g = 255;
	//				}
	//			}

	//			tex[n] = r;
	//			tex[n + 1] = g;
	//			tex[n + 2] = b;
	//			tex[n + 3] = 255;
	//			n += 4;
	//		}
	//	}
	//}

	#pragma omp for
	for (int n = 0; n < texsize * texsize * texsize * 4; n += 4) {


				int j = (n / 4) % texsize;
				int i = floor(((n / 4) % (texsize * texsize)) / texsize);
				int k = floor((n / 4) / (texsize * texsize));

				//This is good for index sample
				float x = (((float)i / texsize) - 0.5f) * -1000.0f; //invert axis
				float y = (((float)j / texsize) - 0.5f) * 1000.0f;
				float z = (((float)k / texsize) - 0.5f) * 1000.0f;

				//this is good for world coordinates sample
				//float x = (((float)i / texsize) - 0.5f) * 100.0f;
				//float y = (((float)j / texsize) - 0.5f) * 100.0f + 25; //250 to put the bunny in the center of the render
				//float z = (((float)k / texsize) - 0.5f) * 100.0f;

				int r = 0;
				int g = 0;
				int b = 0;

				float data = grid->getValue(y - gridy * 1000, x - gridx * 1000, z - gridz * 1000);

				//if (data > 0)
				//std::cout << "value: " << data << "\n";
				if (data != 1.5) //why 1.5?
				{
					r = 255;
					g = 255;
				}

				if (BOUNDINGBOX == 1)
				{
					if ((i == 0 || i == (texsize - 1)) && (k == 0 || k == (texsize - 1)) ||
						(j == 0 || j == (texsize - 1)) && (k == 0 || k == (texsize - 1)) ||
						(i == 0 && j == 0) || (i == 0 && j == (texsize - 1)) ||
						(i == (texsize - 1) && j == 0) || (i == (texsize - 1) && j == (texsize - 1))
						)
					{
						r = 255;
						g = 255;
					}
				}

				tex[n] = r;
				tex[n + 1] = g;
				tex[n + 2] = b;
				tex[n + 3] = 255;
	}

	return tex;
}

void refresh3dTextureAsArray(unsigned char * texture, int gridx, int gridy, int gridz)
{
	int texsize = TEXTURESIZE;

	unsigned char * tex;
	tex = texture;

	//std::cout << " grid x: " << gridx << " y: " << gridy << " z: " << gridz << "  \n";
	//int n = 0;

	//for (int k = 0; k < texsize; k++) {

	//	for (int i = 0; i < texsize; i++)
	//	{
	//		for (int j = 0; j < texsize; j++)
	//		{
	//			//This is good for index sample
	//			float x = (((float)i / texsize) - 0.5f) * -1000.0f; //invert axis
	//			float y = (((float)j / texsize) - 0.5f) * 1000.0f;
	//			float z = (((float)k / texsize) - 0.5f) * 1000.0f;

	//			//this is good for world coordinates sample
	//			//float x = (((float)i / texsize) - 0.5f) * 100.0f;
	//			//float y = (((float)j / texsize) - 0.5f) * 100.0f + 25; //250 to put the bunny in the center of the render
	//			//float z = (((float)k / texsize) - 0.5f) * 100.0f;

	//			int r = 0;
	//			int g = 0;
	//			int b = 0;

	//			//float data = grid->getValue(x - gridx * 1000, y - gridy * 1000, z - gridz * 1000);
	//			float data = grid->getValue(y - gridy * 1000, x - gridx * 1000, z - gridz * 1000);

	//			//if (data > 0)
	//			//if (data != BACKGROUNDVALUE)
	//			if (data != 1.5)
	//			{
	//				r = 255;
	//				g = 255;
	//			}

	//			if (BOUNDINGBOX == 1)
	//			{
	//				if ((i == 0 || i == (texsize - 1)) && (k == 0 || k == (texsize - 1)) ||
	//					(j == 0 || j == (texsize - 1)) && (k == 0 || k == (texsize - 1)) ||
	//					(i == 0 && j == 0) || (i == 0 && j == (texsize - 1)) ||
	//					(i == (texsize - 1) && j == 0) || (i == (texsize - 1) && j == (texsize - 1))
	//					)
	//				{
	//					r = 255;
	//					g = 255;
	//				}
	//			}

	//			tex[n] = r;
	//			tex[n + 1] = g;
	//			tex[n + 2] = b;
	//			tex[n + 3] = 255;
	//			n += 4;
	//		}
	//	}
	//}

	int n = 0;
	
	#pragma omp parallel for shared(tex, grid)
	for (int n = 0; n < texsize * texsize * texsize * 4; n+=4) {

				
				int j = (n / 4) % texsize ;
				int i = floor(((n / 4) % (texsize * texsize))/ texsize);
				int k = floor((n / 4) / (texsize * texsize));
				
				//This is good for index sample
				float x = (((float)i / texsize) - 0.5f) * -1000.0f; //invert axis
				float y = (((float)j / texsize) - 0.5f) * 1000.0f;
				float z = (((float)k / texsize) - 0.5f) * 1000.0f;

				//this is good for world coordinates sample
				//float x = (((float)i / texsize) - 0.5f) * 100.0f;
				//float y = (((float)j / texsize) - 0.5f) * 100.0f + 25; //250 to put the bunny in the center of the render
				//float z = (((float)k / texsize) - 0.5f) * 100.0f;

				int r = 0;
				int g = 0;
				int b = 0;

				//float data = grid->getValue(x - gridx * 1000, y - gridy * 1000, z - gridz * 1000);
				float data = grid->getValue(y - gridy * 1000, x - gridx * 1000, z - gridz * 1000);

				//if (data > 0)
				//if (data != BACKGROUNDVALUE)
				if (data != 1.5)
				{
					r = 255;
					g = 255;
				}

				if (BOUNDINGBOX == 1)
				{
					if ((i == 0 || i == (texsize - 1)) && (k == 0 || k == (texsize - 1)) ||
						(j == 0 || j == (texsize - 1)) && (k == 0 || k == (texsize - 1)) ||
						(i == 0 && j == 0) || (i == 0 && j == (texsize - 1)) ||
						(i == (texsize - 1) && j == 0) || (i == (texsize - 1) && j == (texsize - 1))
						)
					{
						r = 255;
						g = 255;
					}
				}
				#pragma omp critical
				tex[n] = r;
				tex[n + 1] = g;
				tex[n + 2] = b;
				tex[n + 3] = 255;
	}
}

PT(Texture) Render3dBigTexture()
{
	int texsize = TEXTURESIZE;

	texsize *= ((GRIDEXTENSION * 2) + 1);
	PT(Texture)  bunn = new Texture("bunn");
	bunn->setup_3d_texture(texsize, texsize, texsize, Texture::ComponentType::T_float, Texture::Format::F_rgba8);

	//std::cout << " grid x: " << gridx << " y: " << gridy << " z: " << gridz << "  \n";

	for (int k = 0; k < texsize; k++) {
		PNMImage* pPNMImage = new PNMImage(texsize, texsize, 4);

		for (int i = 0; i < texsize; i++)
		{
			for (int j = 0; j < texsize; j++)
			{
				float r = 0.0f;
				float g = 0.0f;
				float b = 0.0f;

				if (BOUNDINGBOX == 1)
				{
					if ((i == 0 || i == (texsize - 1)) && (k == 0 || k == (texsize - 1)) ||
						(j == 0 || j == (texsize - 1)) && (k == 0 || k == (texsize - 1)) ||
						(i == 0 && j == 0) || (i == 0 && j == (texsize - 1)) ||
						(i == (texsize - 1) && j == 0) || (i == (texsize - 1) && j == (texsize - 1))
						)
					{
						r = 1.0f;
						g = 1.0f;
					}
				}

				//int flipy = -j + texsize - 1;
				//std::cout << flipy << "\n";

				pPNMImage->set_red(i, j, r);
				pPNMImage->set_green(i, j, g);
				pPNMImage->set_blue(i, j, b);
				pPNMImage->set_alpha(i, j, 1.0f);
			}
		}

		bunn->load(*pPNMImage, k, 0);

		delete pPNMImage;
	}

	//bunn->write(Filename("woodgrain-#.png"), 0, 0, true, false);

	bunn->set_wrap_u(SamplerState::WrapMode::WM_border_color);
	bunn->set_wrap_v(SamplerState::WrapMode::WM_border_color);
	bunn->set_wrap_w(SamplerState::WrapMode::WM_border_color);
	bunn->set_border_color(LColor(0.0, 0.0, 0.0, 0.0));
	bunn->set_keep_ram_image(true);
	return bunn;
}

void initOffsetVectors()
{
	float j = GRIDEXTENSION * 1.0;
	for (int i = 0; i < ((GRIDEXTENSION * 2) + 1); i++)
	{
		offsetvectorx[i] = j;
		offsetvectory[i] = j;
		offsetvectorz[i] = j;
		j -= 1.0;
	}

	j = -GRIDEXTENSION * 1.0;
	for (int i = 0; i < ((GRIDEXTENSION * 2) + 1); i++)
	{
		offsetvectorxO[i] = j;
		offsetvectoryO[i] = j;
		offsetvectorzO[i] = j;
		j += 1.0;
	}

	int n = -GRIDEXTENSION;

	for (int i = 0; i < ((GRIDEXTENSION * 2) + 1); i++)
	{
		gridoffsetx[i] = n;
		gridoffsety[i] = n;
		gridoffsetz[i] = n;
		n += 1;
	}

}

void refresh3dTexture()
{
	int textsize = TEXTURESIZE * TEXTURESIZE * TEXTURESIZE * 4;

	//PTA_uchar image = bunn->modify_ram_image();
	KeyTriple keycenter = std::make_tuple(0, 0, 0);
	unsigned char * image = gridFrustrum[keycenter];

	int z = 0;
	int y = (TEXTURESIZE - 1);
	int x = 0;

	for (int i = 0; i < textsize; i += 4)
	{
		float x0 = (((float)x / TEXTURESIZE) - 0.5f) * 1000.0f;
		float y0 = (((float)y / TEXTURESIZE) - 0.5f) * 1000.0f; //250 to put the bunny in the center of the render
		float z0 = (((float)z / TEXTURESIZE) - 0.5f) * 1000.0f;

		//rotatey(angledegrees, &x0, &y0, &z0);
		float data = grid->getValue(x0, y0, z0);

		if (data > 0)
		{
			//(uint)image[i] BLUE COLOR
			image[i + 1] = 255; // GREEN COLOR
			image[i + 2] = 255; // RED COLOR
		}
		else
		{
			image[i] = 0;  // BLUE COLOR
			image[i + 1] = 0; // GREEN COLOR
			image[i + 2] = 0; // RED COLOR
		}

		x++;

		if (x == TEXTURESIZE)
		{
			x = 0;
			y--;

			if (y == -1)
			{
				y = (TEXTURESIZE - 1);
				z++;
			}
		}
	}
}

void refresh3dTexture(unsigned char * texture, int gridx, int gridy, int gridz)
{
	int textsize = TEXTURESIZE * TEXTURESIZE * TEXTURESIZE * 4;

	unsigned char * image = texture;

	int z = 0;
	int y = (TEXTURESIZE - 1);
	int x = 0;
	std::cout << "refresh: " << gridx << " " << gridy << " " << gridz << "\n";
	for (int i = 0; i < textsize; i += 4)
	{
		float x0 = (((float)x / TEXTURESIZE) - 0.5f) * 1000.0f;
		float y0 = (((float)y / TEXTURESIZE) - 0.5f) * 1000.0f; //250 to put the bunny in the center of the render
		float z0 = (((float)z / TEXTURESIZE) - 0.5f) * 1000.0f;

		//rotatey(angledegrees, &x0, &y0, &z0);
		//float data = grid->getValue(x0, y0, z0);
		float data = grid->getValue(x0 - gridx * 1000, y0 - gridy * 1000, z0 - gridz * 1000);

		if (data > 0)
		{
			//(uint)image[i] BLUE COLOR
			image[i + 1] = 255; // GREEN COLOR
			image[i + 2] = 255; // RED COLOR
		}
		else
		{
			image[i] = 0;  // BLUE COLOR
			image[i + 1] = 0; // GREEN COLOR
			image[i + 2] = 0; // RED COLOR
		}

		x++;

		if (x == TEXTURESIZE)
		{
			x = 0;
			y--;

			if (y == -1)
			{
				y = (TEXTURESIZE - 1);
				z++;
			}
		}
	}
}

AsyncTask::DoneStatus cameraMotionTask(GenericAsyncTask *task, void *data) {
	// Calculate the new position and orientation (inefficient - change me!)
	double time = globalClock->get_real_time();
	double angledegrees = time * 6.0;
	double angleradians = angledegrees * (3.14 / 180.0);
	//camera.set_pos(1 * sin(angleradians), 0, 1*cos(angleradians) + CAM_z);//orbit around center

	if (FORWARD)
	{
		LVector3f fw = mainWindow->get_render().get_relative_point(camera, LVector3f(0, 0, -1.0));
		fw = fw - camera.get_pos();

		VELX = VELX + ACCELERATION * globalClock->get_dt() * fw.get_x();
		VELY = VELY + ACCELERATION * globalClock->get_dt() * fw.get_z();
		VELZ = VELZ + ACCELERATION * globalClock->get_dt() * fw.get_y();

	}

	if (BACKWARDS)
	{
		LVector3f fw = mainWindow->get_render().get_relative_point(camera, LVector3f(0, 0, 1.0));
		fw = fw - camera.get_pos();

		VELX = VELX + ACCELERATION * globalClock->get_dt() * fw.get_x();
		VELY = VELY + ACCELERATION * globalClock->get_dt() * fw.get_z();
		VELZ = VELZ + ACCELERATION * globalClock->get_dt() * fw.get_y();

	}

	if (LEFT)
	{
		LVector3f fw = mainWindow->get_render().get_relative_point(camera, LVector3f(-1.0, 0, 0.0));
		fw = fw - camera.get_pos();

		VELX = VELX + ACCELERATION * globalClock->get_dt() * fw.get_x();
		VELY = VELY + ACCELERATION * globalClock->get_dt() * fw.get_z();
		VELZ = VELZ + ACCELERATION * globalClock->get_dt() * fw.get_y();
	}

	if (RIGHT)
	{
		LVector3f fw = mainWindow->get_render().get_relative_point(camera, LVector3f(1.0, 0, 0.0));
		fw = fw - camera.get_pos();

		VELX = VELX + ACCELERATION * globalClock->get_dt() * fw.get_x();
		VELY = VELY + ACCELERATION * globalClock->get_dt() * fw.get_z();
		VELZ = VELZ + ACCELERATION * globalClock->get_dt() * fw.get_y();

	}

	if (UP)
	{
		LVector3f fw = mainWindow->get_render().get_relative_point(camera, LVector3f(0.0, -1.0, 0.0));
		fw = fw - camera.get_pos();

		VELX = VELX + ACCELERATION * globalClock->get_dt() * fw.get_x();
		VELY = VELY + ACCELERATION * globalClock->get_dt() * fw.get_z();
		VELZ = VELZ + ACCELERATION * globalClock->get_dt() * fw.get_y();

	}

	if (DOWN)
	{
		LVector3f fw = mainWindow->get_render().get_relative_point(camera, LVector3f(0.0, 1.0, 0.0));
		fw = fw - camera.get_pos();

		VELX = VELX + ACCELERATION * globalClock->get_dt() * fw.get_x();
		VELY = VELY + ACCELERATION * globalClock->get_dt() * fw.get_z();
		VELZ = VELZ + ACCELERATION * globalClock->get_dt() * fw.get_y();

	}

	if (SPINL)
	{
		SPINVEL = SPINVEL + ACCELERATION * globalClock->get_dt() * 50;
	}

	if (SPINR)
	{
		SPINVEL = SPINVEL - ACCELERATION * globalClock->get_dt() * 50;
	}

	VELX -= (VELX  * (1.0 - FRICTION)) * globalClock->get_dt() * 10.0;
	VELY -= (VELY  * (1.0 - FRICTION)) * globalClock->get_dt() * 10.0;
	VELZ -= (VELZ  * (1.0 - FRICTION)) * globalClock->get_dt() * 10.0;
	SPINVEL -= (SPINVEL  * (1.0 - FRICTION)) * globalClock->get_dt() * 10.0;


	if (abs(VELX) < 0.0001)
	{
		VELX = 0.0f;
	}

	if (abs(VELY) < 0.0001)
	{
		VELY = 0.0f;
	}

	if (abs(VELZ) < 0.0001)
	{
		VELZ = 0.0f;
	}

	if (abs(SPINVEL) < 0.0001)
	{
		SPINVEL = 0.0f;
	}


	float nCAM_x = CAM_x - VELX * globalClock->get_dt();
	float nCAM_z = CAM_z - VELY * globalClock->get_dt();
	float nCAM_y = CAM_y - VELZ * globalClock->get_dt();
	ANGLEDEGREES = ANGLEDEGREES + SPINVEL * globalClock->get_dt();


	if (nCAM_x > 0.5 * GRID_SCALE || nCAM_x < -0.5 * GRID_SCALE ||
		nCAM_y > 0.5 * GRID_SCALE || nCAM_y < -0.5 * GRID_SCALE ||
		nCAM_z > 0.5 * GRID_SCALE || nCAM_z < -0.5 * GRID_SCALE
		)
	{
		if (nCAM_x > 0.5 * GRID_SCALE)
		{
			GRID_x--;
			nCAM_x = -0.5 * GRID_SCALE + abs(0.5 * GRID_SCALE - nCAM_x);
		}
		else if (nCAM_x < -0.5 * GRID_SCALE)
		{
			GRID_x++;
			nCAM_x = 0.5 * GRID_SCALE - abs(-0.5 * GRID_SCALE - nCAM_x);
		}

		if (nCAM_y > 0.5 * GRID_SCALE)
		{
			GRID_y--;
			nCAM_y = -0.5 * GRID_SCALE + abs(0.5 * GRID_SCALE - nCAM_y);
		}
		else if (nCAM_y < -0.5 * GRID_SCALE)
		{
			GRID_y++;
			nCAM_y = 0.5 * GRID_SCALE + abs(-0.5 * GRID_SCALE - nCAM_y);
		}

		if (nCAM_z > 0.5 * GRID_SCALE)
		{
			GRID_z--;
			nCAM_z = -0.5 * GRID_SCALE + abs(0.5 * GRID_SCALE - nCAM_z);
		}
		else if (nCAM_z < -0.5 * GRID_SCALE)
		{
			GRID_z++;
			nCAM_z = 0.5 * GRID_SCALE - abs(-0.5 * GRID_SCALE - nCAM_z);
		}

		//Copy/Refresh 3d sectors
		refreshGridFrustrum();
	}

	CAM_x = nCAM_x;
	CAM_z = nCAM_z;
	CAM_y = nCAM_y;

	camera.set_pos(CAM_x, CAM_y, CAM_z);
	camera.set_hpr(0, 0, ANGLEDEGREES);

	LVector3f lookAtDirection = mainWindow->get_render().get_relative_point(camera, LVector3f(0, 0, 1));

	int z = 0;

	//float scalevox = sin(time * 0.05f)*5.0f + 5.0f;
	//std::cout << "sin: " << scalevox << "\n";

	mainQuad[z].set_shader_input("campos", camera.get_pos());
	mainQuad[z].set_shader_input("target", lookAtDirection);
	mainQuad[z].set_shader_input("params", LVector3f(GRID_SCALE, TEXTURE_3D_EXTENSION, INTERNALRES));
	mainQuad[z].set_shader_input("voxparams", LVector3f(VOXEL_SIZE, VOXEL_SIZE, VOXEL_SIZE));
	//mainQuad[z].set_shader_input("params", LVector3f(GRID_SCALE, scalevox, INTERNALRES));

	if (DENOISE == 1)
	{
		mainQuadNorm[z].set_shader_input("campos", camera.get_pos());
		mainQuadNorm[z].set_shader_input("target", lookAtDirection);
		mainQuadNorm[z].set_shader_input("params", LVector3f(GRID_SCALE, TEXTURE_3D_EXTENSION, INTERNALRES));
		mainQuadNorm[z].set_shader_input("voxparams", LVector3f(VOXEL_SIZE, VOXEL_SIZE, VOXEL_SIZE));
		//mainQuadNorm[z].set_shader_input("params", LVector3f(GRID_SCALE, scalevox, INTERNALRES));
	}
	
	//std::cout << "GRID x: " << GRID_x << " y: " << GRID_y << " z: " << GRID_z << "\n";
	//std::cout << "x: " << CAM_x << " y: " << CAM_y << " z: " << CAM_z << "\n";
	// Tell the task manager to continue this task the next frame.


	if (DENOISE == 1)
	{
		//
		////Debug only - noise in texture
		//PTA_uchar image = prevframe_texture->modify_ram_image();

		//for (int yi = 0; yi < INTERNALRES; ++yi) {
		//	for (int xi = 0; xi < INTERNALRES; ++xi) {
		//		int i = (yi * INTERNALRES + xi) * 4;
		//		image[i] = rand() % 256;
		//		image[i + 1] = rand() % 256;
		//		image[i + 2] = rand() % 256;
		//	}
		//}
		////Debug only - noise in texture
		//

		PT(DisplayRegion) displayRegion = mainWindow->get_display_region_3d();

		PT(Texture) tex = displayRegion->get_screenshot();
		CPTA_uchar img = tex->get_ram_image();

		PTA_uchar image = prevframe_texture->modify_ram_image();
		for (int i = 0; i < INTERNALRES * INTERNALRES * 4; i++)
		{
			image[i] = img[i];
		}
		/*for (int yi = 0; yi < INTERNALRES; ++yi) {
			for (int xi = 0; xi < INTERNALRES; ++xi) {
				int i = (yi * INTERNALRES + xi) * 4;
				image[i] = img[i];
				image[i + 1] = img[i + 1];
				image[i + 2] = img[i + 2];
			}
		}*/

		//if (mainWindow->get_graphics_output()) {
		//	mainFramework->get_graphics_engine()->extract_texture_data(prevframe_texture, mainWindow->get_graphics_output()->get_gsg());
			//tempAAbillboard.set_shader_input("data_store", prevframe_texture);
			//tempAAbillboard.set_shader_input("params", LVector3f(GRID_SCALE, GRID_SCALE, INTERNALRES));
		//}
	}

	return AsyncTask::DS_cont;
}

//Refresh texture3d with new values from OpenVDB
AsyncTask::DoneStatus cameraTask(GenericAsyncTask *task, void *data) {
	double time = globalClock->get_real_time();
	double angledegrees = time * 6.0;
	double angleradians = angledegrees * (3.14 / 180.0);


	refresh3dTexture();

	// Tell the task manager to continue this task the next frame.
	return AsyncTask::DS_cont;
}


void RefreshTexture(KeyTriple params)
{
	
	int gridx = std::get<0>(params);
	int gridy = std::get<1>(params);
	int gridz = std::get<2>(params);

	//Need to invert x and y axis
	refresh3dTextureAsArray(gridFrustrum[params], gridy, gridx, gridz);

	renderQueue.push(params);
}


//Dot rasterizer spinning
AsyncTask::DoneStatus spinRasterizerTask(GenericAsyncTask *task, void *data) {
	// Calculate the new position and orientation (inefficient - change me!)
	double time = globalClock->get_real_time();
	double angledegrees = time * 24.0;
	double angleradians = angledegrees * (3.14 / 180.0);

	//const auto startTime = high_resolution_clock::now();

	GeomVertexWriter color(vdata, "color");

#pragma omp parallel for
	for (int i = 0; i < WIDTH * HEIGHT; i++)
	{
		int x = i % WIDTH;
		int y = (int)(i / WIDTH);

		//Inigo quilez camera code
		LVector2f pixel = LVector2f((2.0f * x - WIDTH) / (float)HEIGHT,
			(2.0f * y - HEIGHT) / (float)HEIGHT);

		LVector3f rayorigin = LVector3f(1000.0f * sin(angleradians),
			200.0f,
			1000.0f * cos(angleradians));
		LVector3f target = LVector3f(0, 400.0f, 0);

		//camera matrix
		LVector3f ww = target - rayorigin;
		ww.normalize();

		LVector3f uu = ww.cross(LVector3f(0, 1, 0));
		uu.normalize();

		LVector3f vv = uu.cross(ww);
		vv.normalize();

		//Create view ray
		LVector3f px = LVector3f(pixel.get_x() * uu.get_x(),
			pixel.get_x() * uu.get_y(),
			pixel.get_x() * uu.get_z());
		LVector3f py = LVector3f(pixel.get_y() * vv.get_x(),
			pixel.get_y() * vv.get_y(),
			pixel.get_y() * vv.get_z());
		LVector3f pz = LVector3f(1.5f * ww.get_x(), //1.5f is FOV
			1.5f * ww.get_y(),
			1.5f * ww.get_z());


		LVector3f raydirection = px + py + pz;
		raydirection.normalize();

		//raymarch
		float r = 0.0f;
		float g = 0.0f;
		float b = 0.0f;
		float delta = 0.0f;

		for (int n = 0; n < 25; n++)
		{
			LVector3f pixelsolver = rayorigin + delta * raydirection;

			float data = grid->getValue(pixelsolver.get_x(), pixelsolver.get_y(), pixelsolver.get_z());

			if (data > 0)
			{
				//r = 1.0f;
				//r += data;
				r = 1.0f - 4 * (delta / 5000);
				g = 1.0f - 4 * (delta / 5000);
				break;
			}
			else
			{
				delta = delta + 50.0f;
			}
		}
		color.set_row(i);
		color.set_data4(r, g, b, 1);
	}


	/*
	Normal raycast no parallel execution

	GeomVertexWriter color(vdata, "color");
	for (int i = 0; i < HEIGHT; i++)
	{
		for (int j = 0; j < WIDTH; j++)
		{
			//Inigo quilez camera code
			LVector2f pixel = LVector2f((2.0f * j - WIDTH)  / (float) HEIGHT,
										(2.0f * i - HEIGHT) / (float) HEIGHT);

			LVector3f rayorigin = LVector3f(1000.0f * sin(angleradians),
											200.0f,
											1000.0f * cos(angleradians));
			LVector3f target = LVector3f(0, 400.0f, 0);

			//camera matrix
			LVector3f ww = target - rayorigin;
			ww.normalize();

			LVector3f uu = ww.cross(LVector3f(0, 1, 0));
			uu.normalize();

			LVector3f vv = uu.cross(ww);
			vv.normalize();

			//Create view ray
			LVector3f px = LVector3f(	pixel.get_x() * uu.get_x(),
										pixel.get_x() * uu.get_y(),
										pixel.get_x() * uu.get_z());
			LVector3f py = LVector3f(	pixel.get_y() * vv.get_x(),
										pixel.get_y() * vv.get_y(),
										pixel.get_y() * vv.get_z());
			LVector3f pz = LVector3f(	1.5f * ww.get_x(), //1.5f is FOV
										1.5f * ww.get_y(),
										1.5f * ww.get_z());


			LVector3f raydirection = px + py + pz;
			raydirection.normalize();

			//raymarch
			float r = 0.0f;
			float g = 0.0f;
			float b = 0.0f;
			float delta = 0.0f;

			for (int n = 0; n < 25; n++)
			{
				LVector3f pixelsolver = rayorigin +  delta * raydirection;

				float data = grid->getValue(pixelsolver.get_x(), pixelsolver.get_y(), pixelsolver.get_z());

				if (data > 0)
				{
					//r = 1.0f;
					//r += data;
					r = 1.0f - 4*(delta / 5000);
					g = 1.0f - 4 * (delta / 5000);
					break;
				}
				else
				{
					delta = delta + 50.0f;
				}
			}
			color.set_data4(r, g, b, 1);

		}
	}
	*/
	//const auto endTime = high_resolution_clock::now();

	//print_results("Serial", startTime, endTime);
	// Tell the task manager to continue this task the next frame.
	return AsyncTask::DS_cont;
}

void initGridFrustrum()
{
	KeyTriple key[((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)];

	int z = 0;
	for (int w = 0; w < ((GRIDEXTENSION * 2) + 1); w++)
	{
		for (int v = 0; v < ((GRIDEXTENSION * 2) + 1); v++)
		{
			for (int u = 0; u < ((GRIDEXTENSION * 2) + 1); u++)
			{
				key[z] = std::make_tuple(GRID_y + gridoffsety[v], GRID_x + gridoffsetx[u], GRID_z + gridoffsetz[w]);

				gridTextureArray[z] = Render3dTextureAsArray(gridoffsetx[u], gridoffsety[v], gridoffsetz[w]);
				gridFrustrum[key[z]] = gridTextureArray[z];

				usedTexturesVector[gridTextureArray[z]] = true;

				z++;
			}
		}
	}

	//Only use just 1 big 3d texture
	//KeyTriple centerkey = std::make_tuple(0, 0, 0);
	//Texture * finaltexture = gridFrustrum[centerkey];
	PT(Texture) finaltexture = Render3dBigTexture();

	mainQuad[0].set_texture(finaltexture);

	if (DENOISE == 1) {
		mainQuadNorm[0].set_texture(finaltexture);
	}
}


void initBigTexture()
{
	KeyTriple key[((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)];

	for (int w = 0; w < ((GRIDEXTENSION * 2) + 1); w++)
	{
		for (int v = 0; v < ((GRIDEXTENSION * 2) + 1); v++)
		{
			for (int u = 0; u < ((GRIDEXTENSION * 2) + 1); u++)
			{
				callOpenGLSubImage(GRID_x + gridoffsetx[u], GRID_y + gridoffsety[v], GRID_z + gridoffsetz[w], 0, 0);

				//if (DENOISE == 1) {
				//	callOpenGLSubImage(GRID_x + gridoffsetx[u], GRID_y + gridoffsety[v], GRID_z + gridoffsetz[w], 0, 1);
				//}

			}
		}
	}
}


void refreshGridFrustrum()
{
	GridFrustrum cache = gridFrustrum;
	gridFrustrum.clear();

	KeyTriple key[((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)];

	std::vector<CopyTuple> copyarray;

	//std::vector<RefreshTuple> refresharray;
	std::vector<KeyTriple> refresharray;

	KeyTriple keyold;

	//init used textures vector
	UsedTextures::iterator it;

	for (it = usedTexturesVector.begin(); it != usedTexturesVector.end(); it++)
	{
		usedTexturesVector[it->first] = false;
	}

	int z = 0;
	for (int w = 0; w < ((GRIDEXTENSION * 2) + 1); w++)
	{
		for (int v = 0; v < ((GRIDEXTENSION * 2) + 1); v++)
		{
			for (int u = 0; u < ((GRIDEXTENSION * 2) + 1); u++)
			{
				//invert x and y axis
				key[z] = std::make_tuple(GRID_x + gridoffsety[v], GRID_y + gridoffsetx[u], GRID_z + gridoffsetz[w]); //inverted axis

				if (cache.count(key[z]) == 1)
				{
					gridFrustrum[key[z]] = cache[key[z]];
					usedTexturesVector[cache[key[z]]] = true;
					callOpenGLSubImage(GRID_x + gridoffsety[v], GRID_y + gridoffsetx[u], GRID_z + gridoffsetz[w], 0, 0);

					//if (DENOISE == 1)
					//{
					//	callOpenGLSubImage(GRID_x + gridoffsety[v], GRID_y + gridoffsetx[u], GRID_z + gridoffsetz[w], 0, 1);
					//}
	
				}
				else
				{
					refresharray.push_back(std::make_tuple(GRID_x + gridoffsety[v], GRID_y + gridoffsetx[u], GRID_z + gridoffsetz[w]));
				}

				z++;
			}
		}
	}


	////-----------debugging just refreshing everything----------------------
	//refresharray.push_back(std::make_tuple(0, GRID_x - 1, GRID_y - 1, GRID_z - 1));
	//refresharray.push_back(std::make_tuple(1, GRID_x + 0, GRID_y - 1, GRID_z - 1));
	//refresharray.push_back(std::make_tuple(2, GRID_x + 1, GRID_y - 1, GRID_z - 1));

	//refresharray.push_back(std::make_tuple(3, GRID_x - 1, GRID_y, GRID_z - 1));
	//refresharray.push_back(std::make_tuple(4, GRID_x + 0, GRID_y, GRID_z - 1));
	//refresharray.push_back(std::make_tuple(5, GRID_x + 1, GRID_y, GRID_z - 1));

	//refresharray.push_back(std::make_tuple(6, GRID_x - 1, GRID_y + 1, GRID_z - 1));
	//refresharray.push_back(std::make_tuple(7, GRID_x + 0, GRID_y + 1, GRID_z - 1));
	//refresharray.push_back(std::make_tuple(8, GRID_x + 1, GRID_y + 1, GRID_z - 1));

	//refresharray.push_back(std::make_tuple(9, GRID_x - 1, GRID_y - 1, GRID_z));
	//refresharray.push_back(std::make_tuple(10, GRID_x + 0, GRID_y - 1, GRID_z));
	//refresharray.push_back(std::make_tuple(11, GRID_x + 1, GRID_y - 1, GRID_z));

	//refresharray.push_back(std::make_tuple(12, GRID_x - 1, GRID_y, GRID_z));
	//refresharray.push_back(std::make_tuple(13, GRID_x, GRID_y, GRID_z)); //Dont refresh center, just copy it from near grids
	//refresharray.push_back(std::make_tuple(14, GRID_x + 1, GRID_y, GRID_z));

	//refresharray.push_back(std::make_tuple(15, GRID_x - 1, GRID_y + 1, GRID_z));
	//refresharray.push_back(std::make_tuple(16, GRID_x + 0, GRID_y + 1, GRID_z));
	//refresharray.push_back(std::make_tuple(17, GRID_x + 1, GRID_y + 1, GRID_z));

	//refresharray.push_back(std::make_tuple(18, GRID_x - 1, GRID_y - 1, GRID_z + 1));
	//refresharray.push_back(std::make_tuple(19, GRID_x + 0, GRID_y - 1, GRID_z + 1));
	//refresharray.push_back(std::make_tuple(20, GRID_x + 1, GRID_y - 1, GRID_z + 1));

	//refresharray.push_back(std::make_tuple(21, GRID_x - 1, GRID_y, GRID_z + 1));
	//refresharray.push_back(std::make_tuple(22, GRID_x + 0, GRID_y, GRID_z + 1));
	//refresharray.push_back(std::make_tuple(23, GRID_x + 1, GRID_y, GRID_z + 1));

	//refresharray.push_back(std::make_tuple(24, GRID_x - 1, GRID_y + 1, GRID_z + 1));
	//refresharray.push_back(std::make_tuple(25, GRID_x + 0, GRID_y + 1, GRID_z + 1));
	//refresharray.push_back(std::make_tuple(26, GRID_x + 1, GRID_y + 1, GRID_z + 1));

	for (int i = 0; i < refresharray.size(); i++)
	{
		for (it = usedTexturesVector.begin(); it != usedTexturesVector.end(); it++)
		{
			if (usedTexturesVector[it->first] == false)
			{
				KeyTriple keydata = refresharray[i];
				gridFrustrum[keydata] = it->first;
				usedTexturesVector[it->first] = true;

				refreshQueue.push(keydata);
				callOpenGLSubImage(std::get<0>(keydata), std::get<1>(keydata), std::get<2>(keydata), 1, 0); //refresh texture with empty

				//if (DENOISE == 1)
				//{
				//	callOpenGLSubImage(std::get<0>(keydata), std::get<1>(keydata), std::get<2>(keydata), 1, 1); //refresh texture with empty
				//}

				break;
			}
		}
	}

}

void CleanTexture(PT(Texture) origin)
{
	int textsize = TEXTURESIZE * TEXTURESIZE * TEXTURESIZE * 4;

	PTA_uchar imageorg = origin->modify_ram_image();

	int z = 0;
	int y = (TEXTURESIZE - 1);
	int x = 0;

	for (int i = 0; i < textsize; i += 4)
	{
		imageorg[i] = 0;  // BLUE COLOR
		imageorg[i + 1] = 0; // GREEN COLOR
		imageorg[i + 2] = 0; // RED COLOR
		imageorg[i + 3] = 0; // ALPHA

		x++;

		if (x == TEXTURESIZE)
		{
			x = 0;
			y--;

			if (y == -1)
			{
				y = (TEXTURESIZE - 1);
				z++;
			}
		}
	}

	//origin->reload();
}

unsigned char * CleanTextureArray()
{
	int texsize = TEXTURESIZE;

	unsigned char * tex;
	tex = (unsigned char *)malloc(texsize * texsize * texsize * 4 * (sizeof(char)));

	int n = 0;

	for (int k = 0; k < texsize; k++) {

		for (int i = 0; i < texsize; i++)
		{
			for (int j = 0; j < texsize; j++)
			{
				int r = 0;
				int g = 0;
				int b = 0;

				if (BOUNDINGBOX == 1)
				{
					if ((i == 0 || i == (texsize - 1)) && (k == 0 || k == (texsize - 1)) ||
						(j == 0 || j == (texsize - 1)) && (k == 0 || k == (texsize - 1)) ||
						(i == 0 && j == 0) || (i == 0 && j == (texsize - 1)) ||
						(i == (texsize - 1) && j == 0) || (i == (texsize - 1) && j == (texsize - 1))
						)
					{
						r = 255;
						g = 255;
					}
				}

				tex[n] = r;
				tex[n + 1] = g;
				tex[n + 2] = b;
				tex[n + 3] = 255;
				n += 4;
			}
		}
	}

	return tex;
}

int main(int argc, char *argv[]) {

	grid = new VDBGrid();
	grid->initGrid();
	//grid->movePoints(); //first move then init accessor
	grid->initFastAccessor();

	if (BUNNY == 0) {
		MakeBunny(argc, argv);
	}
	else
	{
		MakeShadertoy(argc, argv);
	}

	return (0);
}

void MakeBunny(int argc, char *argv[])
{
	// Open a new window framework
	PandaFramework framework;
	framework.open_framework(argc, argv);

	load_prc_file_data("", "show-frame-rate-meter 1");

	// Set the window title and open the window
	framework.set_window_title("Dnttg V2 - CEV");
	WindowFramework *window = framework.open_window();
	// Get the camera and store it in a variable.
	camera = window->get_camera_group();
	//window->setup_trackball();//move camera with mouse, errors while trying to move the camera from code directly
	camera.set_pos(0, -3, 0);
	camera.look_at(0, 0, 0);

	//PT(FrameRateMeter) meter = new FrameRateMeter("fps");
	//meter->setup_window(window->get_graphics_output());

	// Load the environment model.
	//NodePath scene = window->load_model(framework.get_models(), "models/environment");
	// Reparent the model to render.
	//scene.reparent_to(window->get_render());
	// Apply scale and position transforms to the model.
	//scene.set_scale(0.15f, 0.15f, 0.15f);
	//scene.set_pos(-8, 42, -1);

	// Load our panda
	//NodePath pandaActor = window->load_model(framework.get_models(), "models/panda-model");
	//pandaActor.set_scale(0.005);
	//pandaActor.reparent_to(window->get_render());

	// Load the walk animation
	//window->load_model(pandaActor, "models/panda-walk4");
	//window->loop_animations(0); // bind models and animations
								//set animations to loop

	/*
	// Create the lerp intervals needed to walk back and forth
	PT(CLerpNodePathInterval) pandaPosInterval1, pandaPosInterval2,
		pandaHprInterval1, pandaHprInterval2;
	pandaPosInterval1 = new CLerpNodePathInterval("pandaPosInterval1",
		13.0, CLerpInterval::BT_no_blend,
		true, false, pandaActor, NodePath());
	pandaPosInterval1->set_start_pos(LPoint3(0, 10, 0));
	pandaPosInterval1->set_end_pos(LPoint3(0, -10, 0));

	pandaPosInterval2 = new CLerpNodePathInterval("pandaPosInterval2",
		13.0, CLerpInterval::BT_no_blend,
		true, false, pandaActor, NodePath());
	pandaPosInterval2->set_start_pos(LPoint3(0, -10, 0));
	pandaPosInterval2->set_end_pos(LPoint3(0, 10, 0));

	pandaHprInterval1 = new CLerpNodePathInterval("pandaHprInterval1", 3.0,
		CLerpInterval::BT_no_blend,
		true, false, pandaActor, NodePath());
	pandaHprInterval1->set_start_hpr(LPoint3(0, 0, 0));
	pandaHprInterval1->set_end_hpr(LPoint3(180, 0, 0));

	pandaHprInterval2 = new CLerpNodePathInterval("pandaHprInterval2", 3.0,
		CLerpInterval::BT_no_blend,
		true, false, pandaActor, NodePath());
	pandaHprInterval2->set_start_hpr(LPoint3(180, 0, 0));
	pandaHprInterval2->set_end_hpr(LPoint3(0, 0, 0));

	// Create and play the sequence that coordinates the intervals
	PT(CMetaInterval) pandaPace;
	pandaPace = new CMetaInterval("pandaPace");
	pandaPace->add_c_interval(pandaPosInterval1, 0,
		CMetaInterval::RS_previous_end);
	pandaPace->add_c_interval(pandaHprInterval1, 0,
		CMetaInterval::RS_previous_end);
	pandaPace->add_c_interval(pandaPosInterval2, 0,
		CMetaInterval::RS_previous_end);
	pandaPace->add_c_interval(pandaHprInterval2, 0,
		CMetaInterval::RS_previous_end);
	pandaPace->loop();
	*/


	//float data = grid->getValue(10.5, -100.2, 50.3);
	//std::cout << "data from grid: " << data;

	//Render3dTexture();

	//Procedurally generate a point in space
	vdata = new GeomVertexData("vertex", GeomVertexFormat::get_v3c4(), Geom::UH_static);


	vdata->set_num_rows(WIDTH * HEIGHT);
	GeomVertexWriter vertex(vdata, "vertex");
	GeomVertexWriter color(vdata, "color");

	for (int i = 0; i < HEIGHT; i++)
	{
		for (int j = 0; j < WIDTH; j++)
		{
			float xpos = (j / (float)WIDTH - 0.5f) * 2;
			float ypos = (i / (float)HEIGHT - 0.5f) * 2;
			vertex.add_data3(xpos, 0.0f, ypos);

			float data = grid->getValue(xpos * 500, 0, ypos * 500);
			if (data == 0) {
				color.add_data4(0, 0, 0, 1);
			}
			else
			{
				color.add_data4(1, 0, 0, 1);
			}
		}
	}


	PT(GeomPoints) prim;
	prim = new GeomPoints(GeomEnums::UH_static);

	for (int i = 0; i < (WIDTH * HEIGHT); i++) {
		prim->add_vertex(i);
	}

	prim->close_primitive();

	PT(Geom) geom;
	geom = new Geom(vdata);
	geom->add_primitive(prim);

	PT(GeomNode) node;
	node = new GeomNode("gnode");
	node->add_geom(geom);

	NodePath nodePath = window->get_render().attach_new_node(node);
	//nodePath.set_render_mode_thickness(10.0f);
	nodePath.set_render_mode_thickness(70.0f);

	SceneGraphAnalyzer sga;
	sga.add_node(window->get_render().node());
	sga.write(std::cerr);


	// Add our task.
	taskMgr->add(new GenericAsyncTask("Spins the Raster", &spinRasterizerTask, nullptr));

	// Do the main loop, equal to run() in python
	//framework.main_loop();

	// This is a simpler way to do stuff every frame,
	// if you're too lazy to create a task.
	Thread *current_thread = Thread::get_current_thread();
	while (framework.do_frame(current_thread)) {
		// Step the interval manager
		CIntervalManager::get_global_ptr()->step();
	}

	framework.close_framework();

	delete grid;
}

void GenerateMainBillboard(int width, int height, WindowFramework * window, PT(Texture) texture)
{
	//Procedurally generate a point in space
	PT(GeomVertexData) vdata = new GeomVertexData("vertex", GeomVertexFormat::get_v3t2(), Geom::UH_static);

	vdata->set_num_rows(4);
	GeomVertexWriter vertex(vdata, "vertex");
	GeomVertexWriter texcoord(vdata, "texcoord");

	vertex.add_data3(0.0f, 0.0f, 0.0f);
	texcoord.add_data2(0.0, 1.0);

	vertex.add_data3(width, 0.0f, 0.0f);
	texcoord.add_data2(1.0, 1.0);

	vertex.add_data3(0.0f, 0.0f, -height);
	texcoord.add_data2(0.0, 0.0);

	vertex.add_data3(width, 0.0f, -height);
	texcoord.add_data2(1.0, 0.0);

	PT(GeomTriangles) prim;
	prim = new GeomTriangles(Geom::UH_static);

	prim->add_vertex(0);
	prim->add_vertex(1);
	prim->add_vertex(2);

	prim->add_vertex(2);
	prim->add_vertex(1);
	prim->add_vertex(3);

	prim->close_primitive();

	PT(Geom) geom;
	geom = new Geom(vdata);
	geom->add_primitive(prim);

	PT(GeomNode) node;
	node = new GeomNode("gnode");
	node->add_geom(geom);

	NodePath nodePath = window->get_pixel_2d().attach_new_node(node);

	nodePath.set_texture(texture);
	nodePath.set_bin("fixed", 1600);

	if (DENOISE == 1)
	{
		PT(Shader) tempAA = Shader::load(Shader::ShaderLanguage::SL_GLSL, "shaders/temporalAA.vert", "shaders/temporalAA.frag");
		nodePath.set_shader(tempAA);

		tempAAbillboard = nodePath;
		tempAAbillboard.set_shader_input("data_store", prevframe_texture);
		tempAAbillboard.set_shader_input("params", LVector3f(GRID_SCALE, GRID_SCALE, INTERNALRES));
	}
}


void InitShader(int index, NodePath nodePath, int type)
{
	PT(TextureStage) ts = new TextureStage("ts");
	ts->set_mode(TextureStage::M_modulate);

	if (type == 0)
	{
		PT(Shader) myShader = Shader::load(Shader::ShaderLanguage::SL_GLSL, "shaders/shader.vert", "shaders/shader.frag");
		//nodePath.set_texture(ts, bunn);
		nodePath.set_shader_input("campos", camera.get_pos());
		nodePath.set_shader_input("params", LVector3f(GRID_SCALE, GRID_SCALE, INTERNALRES));
		nodePath.set_shader_input("target", mainWindow->get_render().get_relative_point(camera, LVector3f(0, 0, 1)));
		nodePath.set_shader(myShader);
		nodePath.set_transparency(TransparencyAttrib::Mode::M_alpha);

		mainQuad[index] = nodePath;
	}
	else
	{
		PT(Shader) myShader = Shader::load(Shader::ShaderLanguage::SL_GLSL, "shaders/shadernorm.vert", "shaders/shadernorm.frag");
		//nodePath.set_texture(ts, bunn);
		nodePath.set_shader_input("campos", camera.get_pos());
		nodePath.set_shader_input("params", LVector3f(GRID_SCALE, GRID_SCALE, INTERNALRES));
		nodePath.set_shader_input("target", mainWindow->get_render().get_relative_point(camera, LVector3f(0, 0, 1)));
		nodePath.set_shader(myShader);
		nodePath.set_transparency(TransparencyAttrib::Mode::M_alpha);

		mainQuadNorm[index] = nodePath;
	}

	int u = (int)(index % ((GRIDEXTENSION * 2) + 1)) - GRIDEXTENSION;
	int v = (int)(index / ((GRIDEXTENSION * 2) + 1)) % ((GRIDEXTENSION * 2) + 1) - GRIDEXTENSION;
	int w = (int)(index / (((GRIDEXTENSION * 2) + 1) * ((GRIDEXTENSION * 2) + 1))) - GRIDEXTENSION;

	u = u - GRID_x;
	v = v - GRID_y;
	w = w - GRID_z;

	float distance = sqrtf(u * u + v * v + w * w) * -100.0f + 200.0;

	//std::cout << " z " << index << " distance: " << distance << "\n";

	nodePath.set_bin("fixed", (int)distance);
	//nodePath.set_bin("fixed", (int) (index / ((GRIDEXTENSION * 2) + 1))); //I put 100 in the center quad WARN


}

void GenerateBillboard(int width, int height, WindowFramework * window, int index, bool useBuffer, NodePath parentNode, int centerx, int centery, int type)
{
	//Procedurally generate a point in space
	PT(GeomVertexData) vdatab = new GeomVertexData("vertex", GeomVertexFormat::get_v3c4(), Geom::UH_static);

	vdatab->set_num_rows(4);
	GeomVertexWriter vertex(vdatab, "vertex");
	GeomVertexWriter color(vdatab, "color");

	vertex.add_data3(0.0f, 0.0f, 0.0f);
	color.add_data4(1.0, 0, 0, 1);

	vertex.add_data3(width, 0.0f, 0.0f);
	color.add_data4(0.0, 1.0, 0, 1);

	vertex.add_data3(0.0f, 0.0f, -height);
	color.add_data4(0.0, 0, 1.0f, 1);

	vertex.add_data3(width, 0.0f, -height);
	color.add_data4(1.0, 1.0, 0, 1);

	PT(GeomTriangles) prim;
	prim = new GeomTriangles(Geom::UH_static);

	//order of triangles should be opposite
	prim->add_vertex(2);
	prim->add_vertex(1);
	prim->add_vertex(0);

	prim->add_vertex(3);
	prim->add_vertex(1);
	prim->add_vertex(2);

	prim->close_primitive();

	PT(Geom) geom;
	geom = new Geom(vdatab);
	geom->add_primitive(prim);

	PT(GeomNode) node;
	node = new GeomNode("gnode" + index);
	node->add_geom(geom);

	if (useBuffer)
	{
		NodePath nodePath = parentNode.attach_new_node(node);
		InitShader(index, nodePath, type);

		nodePath.set_pos(centerx, 0, centery);

		if (index > 20)
		{
			nodePath.hide();
		}
	}
	else
	{
		NodePath nodePath = window->get_pixel_2d().attach_new_node(node);
		InitShader(index, nodePath, type);
	}

}

void GeneratePrePassBillboard(int width, int height, WindowFramework * window, NodePath parentNode, int centerx, int centery)
{
	//Procedurally generate a point in space
	PT(GeomVertexData) vdatab = new GeomVertexData("vertex", GeomVertexFormat::get_v3t2(), Geom::UH_static);

	vdatab->set_num_rows(4);
	GeomVertexWriter vertex(vdatab, "vertex");
	GeomVertexWriter texcoord(vdatab, "texcoord");

	vertex.add_data3(0.0f, 0.0f, 0.0f);
	texcoord.add_data2(0.0, 1.0);

	vertex.add_data3(width, 0.0f, 0.0f);
	texcoord.add_data2(1.0, 1.0);

	vertex.add_data3(0.0f, 0.0f, -height);
	texcoord.add_data2(0.0, 0.0);

	vertex.add_data3(width, 0.0f, -height);
	texcoord.add_data2(1.0, 0.0);

	PT(GeomTriangles) prim;
	prim = new GeomTriangles(Geom::UH_static);

	//order of triangles should be opposite
	prim->add_vertex(2);
	prim->add_vertex(1);
	prim->add_vertex(0);

	prim->add_vertex(3);
	prim->add_vertex(1);
	prim->add_vertex(2);

	prim->close_primitive();

	PT(Geom) geom;
	geom = new Geom(vdatab);
	geom->add_primitive(prim);

	PT(GeomNode) node;
	node = new GeomNode("gnodeprepass");
	node->add_geom(geom);

	NodePath nodePath = parentNode.attach_new_node(node);

	nodePath.set_texture(rendertexture01);//diffuse

	PT(Shader) denoise = Shader::load(Shader::ShaderLanguage::SL_GLSL, "shaders/denoise.vert", "shaders/denoise.frag");
	nodePath.set_shader(denoise);
	nodePath.set_transparency(TransparencyAttrib::Mode::M_alpha);
	nodePath.set_shader_input("params", LVector3f(GRID_SCALE, GRID_SCALE, INTERNALRES));

	PT(TextureStage) ts = new TextureStage("ts");
	ts->set_mode(TextureStage::M_modulate);

	nodePath.set_texture(ts, rendertexture02);//normal

	nodePath.set_bin("fixed", 1500); //I put 150 in the center quad WARN, this is render order

	nodePath.set_pos(centerx, 0, centery);
}


void GenerateTextureBuffer(int width, int height, WindowFramework * window, NodePath testscene)
{
	PT(GraphicsOutput) mybuffer;
	PT(Texture) mytexture;
	PT(Camera) mycamera;
	PT(DisplayRegion) region;
	NodePath mycameraNP;
	NodePath myscene;

	int texsize = INTERNALRES;

	mybuffer = window->get_graphics_output()->make_texture_buffer("Main Buffer", texsize, texsize);
	mytexture = mybuffer->get_texture();
	mybuffer->set_sort(-100);
	mycamera = new Camera("main camera");
	mycameraNP = window->get_render().attach_new_node(mycamera);
	region = mybuffer->make_display_region();
	region->set_camera(mycameraNP);
	myscene = NodePath("Main Scene");
	myscene.set_depth_test(false);
	myscene.set_depth_write(false);
	mycameraNP.reparent_to(myscene);

	PT(OrthographicLens) lens = new OrthographicLens();
	lens->set_film_size(width, height);
	lens->set_near_far(-1000, 1000);
	mycamera->set_lens(lens);

	//Debug only - check if texture can be displayer on main billboard
	//LoaderOptions options;
	//PT(Texture) tex = TexturePool::load_texture("p.png", 0, false, options);
	//GenerateMainBillboard(width, height, window, tex);

	//Debug only - test if we can add nodes to myscene
	//testscene.reparent_to(myscene);

	PT(GraphicsOutput) mybuffernorm;
	PT(Texture) mytexturenorm;
	PT(Camera) mycameranorm;
	PT(DisplayRegion) regionnorm;
	NodePath mycameraNPnorm;
	NodePath myscenenorm;

	mybuffernorm = window->get_graphics_output()->make_texture_buffer("Normal Buffer", texsize, texsize);
	mytexturenorm = mybuffernorm->get_texture();
	mybuffernorm->set_sort(-100);
	mycameranorm = new Camera("main camera");
	mycameraNPnorm = window->get_render().attach_new_node(mycameranorm);
	regionnorm = mybuffernorm->make_display_region();
	regionnorm->set_camera(mycameraNPnorm);
	myscenenorm = NodePath("Main Scene");
	myscenenorm.set_depth_test(false);
	myscenenorm.set_depth_write(false);
	mycameraNPnorm.reparent_to(myscenenorm);

	PT(OrthographicLens) lensnorm = new OrthographicLens();
	lensnorm->set_film_size(width, height);
	lensnorm->set_near_far(-1000, 1000);
	mycameranorm->set_lens(lensnorm);


	PT(GraphicsOutput) mybufferaa;
	PT(Texture) mytextureaa;
	PT(Camera) mycameraaa;
	PT(DisplayRegion) regionaa;
	NodePath mycameraNPaa;
	NodePath mysceneaa;

	mybufferaa = window->get_graphics_output()->make_texture_buffer("AA Buffer", texsize, texsize);
	mytextureaa = mybufferaa->get_texture();
	mybufferaa->set_sort(-100);
	mycameraaa = new Camera("main camera");
	mycameraNPaa = window->get_render().attach_new_node(mycameraaa);
	regionaa = mybufferaa->make_display_region();
	regionaa->set_camera(mycameraNPaa);
	mysceneaa = NodePath("Main Scene");
	mysceneaa.set_depth_test(false);
	mysceneaa.set_depth_write(false);
	mycameraNPaa.reparent_to(mysceneaa);

	PT(OrthographicLens) lensaa = new OrthographicLens();
	lensaa->set_film_size(width, height);
	lensaa->set_near_far(-1000, 1000);
	mycameraaa->set_lens(lensaa);

	rendertexture01 = mytexture;
	rendertexture02 = mytexturenorm;
	rendertexture03 = mytextureaa;


	if (DENOISE == 0) {
		GenerateMainBillboard(width, height, window, mytexture);
	}
	else
	{
		GenerateMainBillboard(width, height, window, mytextureaa);
	}

	GenerateBillboard(width, height, window, 0, true, myscene, -width / 2, height / 2, 0);

	if (DENOISE == 1) {
		GenerateBillboard(width, height, window, 0, true, myscenenorm, -width / 2, height / 2, 1);
	}

	if (DENOISE == 1)
	{
		GeneratePrePassBillboard(width, height, window, mysceneaa, -width / 2, height / 2);
	}

}

void MakeShadertoy(int argc, char *argv[])
{
	// Open a new window framework
	PandaFramework framework;
	framework.open_framework(argc, argv);

	load_prc_file_data("", "show-frame-rate-meter 1");

	// Set the window title and open the window
	framework.set_window_title("Dnttg V2 - CEV");
	WindowFramework *window = framework.open_window();

	mainWindow = window;
	mainFramework = &framework;

	WindowProperties wp;
	framework.get_default_window_props(wp);
	int width = wp.get_x_size();
	int height = wp.get_y_size();

	window->enable_keyboard();
	window->get_render().set_antialias(AntialiasAttrib::M_none);

	framework.define_key("arrow_up", "advance cam", advanceCamera, nullptr);
	framework.define_key("arrow_up-up", "advance cam", advanceCameraUp, nullptr);
	framework.define_key("arrow_down", "back cam", backCamera, nullptr);
	framework.define_key("arrow_down-up", "back cam", backCameraUp, nullptr);
	framework.define_key("arrow_left", "sl cam", spinLCamera, nullptr);
	framework.define_key("arrow_left-up", "sl cam", spinLCameraUp, nullptr);
	framework.define_key("arrow_right", "sr cam", spinRCamera, nullptr);
	framework.define_key("arrow_right-up", "sr cam", spinRCameraUp, nullptr);

	framework.define_key("w", "fow cam", advanceCamera, nullptr);
	framework.define_key("w-up", "fow cam", advanceCameraUp, nullptr);
	framework.define_key("s", "back cam", backCamera, nullptr);
	framework.define_key("s-up", "back cam", backCameraUp, nullptr);
	framework.define_key("a", "left cam", leftCamera, nullptr);
	framework.define_key("a-up", "left cam", leftCameraUp, nullptr);
	framework.define_key("d", "right cam", rightCamera, nullptr);
	framework.define_key("d-up", "right cam", rightCameraUp, nullptr);

	framework.define_key("q", "up cam", upCamera, nullptr);
	framework.define_key("q-up", "up cam", upCameraUp, nullptr);

	framework.define_key("e", "down cam", downCamera, nullptr);
	framework.define_key("e-up", "down cam", downCameraUp, nullptr);

	framework.define_key("space", "modify grid", spawnSphere, nullptr);

	if (DENOISE == 1) {
		prevframe_texture = new Texture();
		prevframe_texture->setup_2d_texture(INTERNALRES, INTERNALRES, Texture::T_unsigned_byte, Texture::F_rgba8);
		prevframe_texture->set_clear_color(LColor(0, 0, 0, 0));
	}

	// Get the camera and store it in a variable.
	camera = window->get_camera_group();
	camera.set_pos(0, 0, 1);
	//window->setup_trackball(); //move camera with mouse, errors while trying to move the camera from code directly

	initOffsetVectors();

	bunn = Render3dTexture(0, 0, 0);//dummy bunny for testing
	emptyTexture = Render3dTexture(0, 0, 0);
	emptyTextureArray = CleanTextureArray(); //generate an empty texture array;
	//emptyTextureArray = Render3dTextureAsArray(0, 0, 0); //generate a bunny texture array; just for debugging
	CleanTexture(emptyTexture); //generate an empty texture;


	std::cout << "max textures: " << window->get_graphics_output()->get_gsg()->get_max_texture_stages() << "\n";

	// debug only - Load the environment model and test it on the billboard
	NodePath envscene = window->load_model(framework.get_models(), "models/environment");
	// Apply scale and position transforms to the model.
	envscene.set_scale(0.25f, 0.25f, 0.25f);
	envscene.set_pos(-8, 42, 0);

	//generate compute shader
	//copyTextureShader = Shader::load_compute(Shader::ShaderLanguage::SL_GLSL, "shaders/copytexture.glsl");

	GenerateTextureBuffer(width, height, window, envscene);

	initGridFrustrum();

	SceneGraphAnalyzer sga;
	sga.add_node(window->get_pixel_2d().node());
	sga.write(std::cerr);

	// Add our task.
	// If we specify custom data instead of NULL, it will be passed as the second argument
	// to the task function.
	AsyncTaskChain *chain = taskMgr->make_task_chain("changevdbgrid");
	chain->set_num_threads(1);
	chain->set_thread_priority(ThreadPriority::TP_urgent);

	//renderChain = taskMgr->make_task_chain("renderchain");
	//renderChain->set_num_threads(2);
	//renderChain->set_thread_priority(ThreadPriority::TP_urgent);

	GenericAsyncTask* modifytask = new GenericAsyncTask("Modify Grid", &modifyGrid, nullptr); //modify after spawn
	modifytask->set_task_chain("changevdbgrid");

	GenericAsyncTask* refreshtask = new GenericAsyncTask("Refresh Grid", &refreshGrid, nullptr);
	refreshtask->set_task_chain("changevdbgrid");


	taskMgr->add(new GenericAsyncTask("Camera Motion", &cameraMotionTask, nullptr));
	taskMgr->add(modifytask);
	taskMgr->add(refreshtask);

	int state = 0;
	// This is a simpler way to do stuff every frame,
	// if you're too lazy to create a task.
	Thread *current_thread = Thread::get_current_thread();
	while (framework.do_frame(current_thread)) {
		// Step the interval manager
		CIntervalManager::get_global_ptr()->step();

		//horrible way to load a texture on the first frame
		switch (state)
		{
		case 0:
		{
			if (glTextureSubImage3D == 0) {
				glTextureSubImage3D = (PFNGLTEXTURESUBIMAGE3DPROC)GetAnyGLFuncAddress("glTextureSubImage3D");
			}
			else {
				state = 1;
			}
			break;
		}
		case 1:
		{
			initBigTexture();
			state = 2;
			break;
		}
		default:
		{
			break;
		}


		}

		while (renderQueue.size() > 0) {
			//std::cout << "Queue size: " << renderQueue.size() << "\n";
			KeyTriple args = renderQueue.front();
			renderQueue.pop();

			callOpenGLSubImage(std::get<0>(args), std::get<1>(args), std::get<2>(args), 0, 0); //put last parameter in 0 to render correctly

			//if (DENOISE == 1)
			//{
			//	callOpenGLSubImage(std::get<0>(args), std::get<1>(args), std::get<2>(args), 0, 1); //put last parameter in 0 to render correctly
			//}
		}

	}

	framework.close_framework();

	delete grid;
}

void callOpenGLSubImage(int posx, int posy, int posz, int debug, int quad)
{
	//std::cout << "Texture gl function address: " << (int)glTextureSubImage3D << "\n";
	PT(GraphicsStateGuardianBase) gsg = mainWindow->get_graphics_window()->get_gsg();
	//GLTextureContext *GLtex = DCAST(GLTextureContext, mainQuad[0].get_texture()->prepare_now(0, gsg->get_prepared_objects(), gsg));
	GLTextureContext *GLtex;
	
	if (quad == 0) {
		GLtex = DCAST(GLTextureContext, mainQuad[0].get_texture()->prepare_now(0, gsg->get_prepared_objects(), gsg));
	}
	else
	{
		GLtex = DCAST(GLTextureContext, mainQuadNorm[0].get_texture()->prepare_now(0, gsg->get_prepared_objects(), gsg));
	}

	GLuint texA = GLtex->_index;
	GLuint internal_format = GLtex->_internal_format;
	GLuint wi = GLtex->_width;
	GLuint he = GLtex->_height;
	GLuint de = GLtex->_depth;
	//std::cout << "Texture gl id: " << texA << " int format: " << internal_format << " " << wi << " " << he << " " << de << "\n";


	unsigned char * tex;
	int bigsize = TEXTURESIZE * ((GRIDEXTENSION * 2) + 1);

	int centerx = GRIDEXTENSION;
	int centery = GRIDEXTENSION;
	int centerz = GRIDEXTENSION;


	int finalposx;
	int finalposy;
	int finalposz;


	finalposx = (centerx - posx + GRID_x) * TEXTURESIZE;
	finalposy = (centery + posy - GRID_y) * TEXTURESIZE;
	finalposz = (centerz - posz + GRID_z) * TEXTURESIZE;
	

	KeyTriple key = std::make_tuple(posx, posy, posz);

	if (debug == 0) {
		tex = gridFrustrum[key];
	}
	else
	{
		tex = emptyTextureArray;
	}

	//std::cout << "memory to access: " << (int)tex << "\n";
	//std::cout << " - " << centerx << " " << centery << " " <<  centerz << " " <<  GRID_x << " " <<  GRID_y << " " <<  GRID_z << "\n";
	//std::cout << " x y z " << posx << " " << posy << " " << posz << " " << finalposx << " " << finalposy << " " << finalposz << " " << bigsize << "\n";
	//std::cout << "NANA " << finalposx << " " << finalposy << " " << finalposz << "\n";

	glBindTexture(GL_TEXTURE_3D, texA);
	glTextureSubImage3D(texA, 0, finalposx, finalposy, finalposz, TEXTURESIZE, TEXTURESIZE, TEXTURESIZE, GL_RGBA, GL_UNSIGNED_BYTE, (GLvoid*)tex);

	//GLenum err;
	//err = glGetError();

	//std::cout << "Texture gl error: " << err << "\n";
	//free(tex);
}
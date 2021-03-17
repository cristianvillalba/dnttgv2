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

#include "Dnntgv2.h"


#define WIDTH 200  //for bunny - for first raycaster
#define HEIGHT 200 //for bunny - for first raycaster
#define INTERNALRES 128 //internal texture resolution
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

NodePath mainQuad[((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)];
PT(Texture) gridTextureArray[((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)];
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

//Global camera pos;
float CAM_x = 0.0;
float CAM_y = 0.0;
float CAM_z = 0.0;

//Global camera pos;
int GRID_x = 0;
int GRID_y = 0;
int GRID_z = 0;

//Grid Scale
float GRID_SCALE = 10.0f;

//grid frustrum
GridFrustrum gridFrustrum;

//Velocity vector/Acceleration
float VELX = 0.0;
float VELY = 0.0;
float VELZ = 0.0;
float ACCELERATION = 20.0;
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

//---Support structures
struct TaskArgs{
	int index;
	int gridx;
	int gridy;
	int gridz;
};

struct CopyArgs {
	int index;
	PT(Texture) data;
};


TaskArgs* _alltasks[((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)];
CopyArgs* _copytasks[((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)];

//---Support structures


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
		//LVector3f fw = mainWindow->get_render().get_relative_point(camera, LVector3f(0, 0, -0.10));
		LVector3f fw = mainWindow->get_render().get_relative_point(camera, LVector3f(0, 0, 0));
		fw /= GRID_SCALE;
		fw *= -1000;
		fw.add_x(GRID_x * 1000);
		fw.add_y(GRID_y * 1000);
		fw.add_z(GRID_z * 1000);

		std::cout << "spawn at x: " << fw.get_x() << "spawn at y: " << fw.get_y() << "spawn at z: " << fw.get_z() << "\n";

		fw *= -1; //invert coordinates to spawn
		grid->spawnSphere(fw);

		//refresh3dTexture();
		
		refresh3dTexture(gridTextureArray[13], GRID_x, GRID_y, GRID_z);//only refresh center
		//gridTextureArray[13] = Render3dTexture(GRID_x, GRID_y, GRID_z);

		//KeyTriple key = std::make_tuple(GRID_x, GRID_y, GRID_z);
		//gridFrustrum[key] = gridTextureArray[13];

		//refreshGridFrustrum();

		SPAWN = false;
	}

	return AsyncTask::DS_cont;
}

AsyncTask::DoneStatus refreshGrid(GenericAsyncTask *task, void *data)
{
	if (*pREFRESHGRID)
	{
		mutex.acquire();
		refreshGridFrustrum();
		*pREFRESHGRID = false;
		mutex.release();
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
				float data = grid->getValue(x - gridx * 1000,  y - gridy * 1000, z - gridz * 1000);
	
				if (data > 0 )
				{
					r = 1.0f ;
					g = 1.0f ;
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
	PTA_uchar image = gridFrustrum[keycenter]->modify_ram_image();

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

PT(Texture) RenderShadows(int gridx, int gridy, int gridz)
{
	int texsize = TEXTURESIZE * 2;

	PT(Texture)  bunn = new Texture("bunn");
	bunn->setup_3d_texture(texsize, texsize, texsize, Texture::ComponentType::T_float, Texture::Format::F_rgba8);


	for (int k = 0; k < texsize; k++) {
		PNMImage* pPNMImage = new PNMImage(texsize, texsize, 4);

		for (int i = 0; i < texsize; i++)
		{
			for (int j = 0; j < texsize; j++)
			{
				//This is good for index sample
				float x = (((float)i / texsize) - 0.5f) * 1000.0f * 3.0; //to expand 3 grids across axis
				float y = (((float)j / texsize) - 0.5f) * 1000.0f * 3.0; //to expand 3 grids across axis
				float z = (((float)k / texsize) - 0.5f) * 1000.0f * 3.0; //to expand 3 grids across axis

				float r = 0.0f;
				float g = 0.0f;
				float b = 0.0f;

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

				pPNMImage->set_red(i, j, r);
				pPNMImage->set_green(i, j, g);
				pPNMImage->set_blue(i, j, b);
				pPNMImage->set_alpha(i, j, 1.0f);
			}
		}

		bunn->load(*pPNMImage, k, 0);

		delete pPNMImage;
	}

	//bunn->write(Filename("shadow-#.png"), 0, 0, true, false);

	bunn->set_wrap_u(SamplerState::WrapMode::WM_border_color);
	bunn->set_wrap_v(SamplerState::WrapMode::WM_border_color);
	bunn->set_wrap_w(SamplerState::WrapMode::WM_border_color);
	bunn->set_border_color(LColor(0.0, 0.0, 0.0, 0.0));
	bunn->set_keep_ram_image(true);
	return bunn;
}

void refresh3dTexture(PT(Texture) texture, int gridx, int gridy, int gridz)
{
	int textsize = TEXTURESIZE * TEXTURESIZE * TEXTURESIZE * 4;

	PTA_uchar image = texture->modify_ram_image();

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

	VELX = VELX * 0.9;
	VELY = VELY * 0.9;
	VELZ = VELZ * 0.9;

	SPINVEL = SPINVEL * 0.9;

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

		if (!*pREFRESHGRID)
		{
			mutex.acquire();
			*pREFRESHGRID = true;
			mutex.release();
		}
	}

	//CAM_x = nCAM_x  - floor(nCAM_x + 0.5); //keep cam always in middle range -0.5 0.5
	//CAM_z = nCAM_z  - floor(nCAM_z + 0.5); //keep cam always in middle range -0.5 0.5
	CAM_x = nCAM_x;
	CAM_z = nCAM_z;
	CAM_y = nCAM_y;

	camera.set_pos(CAM_x, CAM_y, CAM_z);
	camera.set_hpr(0, 0, ANGLEDEGREES);

	ZOrdering();

	LVector3f lookAtDirection = mainWindow->get_render().get_relative_point(camera, LVector3f(0, 0, 1));

	//update all quads
	
	int z = 0;
	for (int w = 0; w < ((GRIDEXTENSION * 2) + 1); w++)
	{
		for (int v = 0; v < ((GRIDEXTENSION * 2) + 1); v++)
		{
			for (int u = 0; u < ((GRIDEXTENSION * 2) + 1); u++)
			{
				mainQuad[z].set_shader_input("campos", camera.get_pos() - LVector3f(offsetvectorx[u], offsetvectory[v], offsetvectorz[w]) * GRID_SCALE);
				mainQuad[z].set_shader_input("target", lookAtDirection - LVector3f(offsetvectorx[u], offsetvectory[v], offsetvectorz[w]) * GRID_SCALE);
				mainQuad[z].set_shader_input("params", LVector3f(GRID_SCALE, GRID_SCALE, INTERNALRES));
				z++;
			}
		}
	}
	
		

	
	/*
	if (-floor(CAM_x + 0.5) != GRID_x ||
		-floor(CAM_y + 0.5) != GRID_y ||
		-floor(CAM_z + 0.5) != GRID_z)
	{
		
		GRID_x = -floor(CAM_x + 0.5);
		GRID_y = -floor(CAM_y + 0.5);
		GRID_z = -floor(CAM_z + 0.5);
		
		if (!*pREFRESHGRID)
		{
			mutex.acquire();
			*pREFRESHGRID = true;
			mutex.release();
		}
	}
	*/

	//std::cout << "GRID x: " << GRID_x << " y: " << GRID_y << " z: " << GRID_z << "\n";
	//std::cout << "x: " << CAM_x << " y: " << CAM_y << " z: " << CAM_z << "\n";
	// Tell the task manager to continue this task the next frame.
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

AsyncTask::DoneStatus copyParallelTextureTask(GenericAsyncTask *task, void *data)
{
	CopyArgs* taskdata = (CopyArgs*)data;

	//CopyTexture(taskdata->data, gridTextureArray[taskdata->index]);

	return AsyncTask::DS_done;
}

AsyncTask::DoneStatus renderParallelTextureTask(GenericAsyncTask *task, void *data) {

	TaskArgs* taskdata = (TaskArgs*)data;

	gridTextureArray[taskdata->index] = Render3dTexture(taskdata->gridx, taskdata->gridy, taskdata->gridz);

	//std::cout << "vector index: " << taskdata->index <<"\n";

	KeyTriple key = std::make_tuple(taskdata->gridx, taskdata->gridy, taskdata->gridz);
	gridFrustrum[key] = gridTextureArray[taskdata->index];
	mainQuad[taskdata->index].set_texture(gridFrustrum[key]);

	return AsyncTask::DS_done;
}

bool CompareZPos(Zorder n1, Zorder n2)
{
	LVector3f n1distance = LVector3f(CAM_x, CAM_y, CAM_z) - n1.pos;
	LVector3f n2distance = LVector3f(CAM_x, CAM_y, CAM_z) - n2.pos;
	bool log = false;

	if (n1distance.length() < n2distance.length()) {
		return true;
	}
	else
	{
		return false;
	}
}

void ZOrdering()
{
	std::vector<Zorder> zorder;
	
	Zorder data[((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)];
	LVector3f key[((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)];

	int z = 0;

	for (int w = 0; w < ((GRIDEXTENSION * 2) + 1); w++)
	{
		for (int v = 0; v < ((GRIDEXTENSION * 2) + 1); v++)
		{
			for (int u = 0; u < ((GRIDEXTENSION * 2) + 1); u++)
			{
				key[z] = LVector3f(offsetvectorxO[u]*1000, offsetvectoryO[v] * 1000, offsetvectorzO[w] * 1000);

				data[z].pos = key[z];
				data[z].quad = mainQuad[z];
				data[z].index = z;

				z++;
			}
		}
	}
	
	
	for (int z = 0; z < ((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1); z++)
	{
		zorder.push_back(data[z]);
	}

	std::sort(zorder.begin(), zorder.end(), CompareZPos);

	int i = 0;

	for (std::vector<Zorder>::iterator it = zorder.begin(); it != zorder.end(); ++it) {
		(*it).quad.set_bin("fixed", -i);
		i++;
	}
		
}

void CopyAndRefreshTexture(CopyTuple params, GridFrustrum cache)
{
	//// Create a dummy node and apply the shader to it
	//NodePath dummy("dummy");
	//dummy.set_shader(copyTextureShader);
	//dummy.set_shader_input("fromTex", cache[std::get<1>(params)]);
	//dummy.set_shader_input("toTex", gridTextureArray[std::get<0>(params)]);

	//// Retrieve the underlying ShaderAttrib
	//CPT(ShaderAttrib) sattr = DCAST(ShaderAttrib,
	//	dummy.get_attrib(ShaderAttrib::get_class_type()));

	//// Our image has 32x32 tiles
	//LVecBase3i work_groups(512 / 16, 512 / 16, 1);

	//// Dispatch the compute shader, right now!
	//GraphicsEngine *engine = GraphicsEngine::get_global_ptr();
	//engine->dispatch_compute(work_groups, sattr, mainWindow->get_graphics_window()->get_gsg());
	
	//std::cout << "copy from: " << std::get<0>(std::get<1>(params)) << " " << std::get<1>(std::get<1>(params)) << " " << std::get<2>(std::get<1>(params)) << " to index: " << std::get<0>(params)  << "\n";

	gridTextureArray[std::get<0>(params)] = cache[std::get<1>(params)];

	KeyTriple key = std::get<1>(params);
	gridFrustrum[key] = gridTextureArray[std::get<0>(params)];
	mainQuad[std::get<0>(params)].set_texture(gridTextureArray[std::get<0>(params)], 1);

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
				key[z] = std::make_tuple(GRID_x + gridoffsetx[u], GRID_y + gridoffsety[v], GRID_z + gridoffsetz[w]);

				gridFrustrum[key[z]] = Render3dTexture(gridoffsetx[u], gridoffsety[v], gridoffsetz[w]);

				gridTextureArray[z] = gridFrustrum[key[z]];
				mainQuad[z].set_texture(gridFrustrum[key[z]]);

				z++;
			}
		}
	}

	//array to use them as parametes in parallel tasks
	for (int i = 0; i < ((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1); i++)
	{
		TaskArgs * _taskArgs = new TaskArgs();

		_alltasks[i] = _taskArgs;
	}

	//array to use them as parametes in parallel tasks
	for (int i = 0; i < ((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1); i++)
	{
		CopyArgs * _copyArgs = new CopyArgs();

		_copytasks[i] = _copyArgs;
	}
}

//void refreshGridFrustrum()
//{
//	CopyTexture(bunn, gridTextureArray[3]);
//	mainQuad04.set_texture(gridTextureArray[3]);
//}

void refreshGridFrustrum()
{
	KeyTriple key[((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)];

	int z = 0;
	for (int w = 0; w < ((GRIDEXTENSION * 2) + 1); w++)
	{
		for (int v = 0; v < ((GRIDEXTENSION * 2) + 1); v++)
		{
			for (int u = 0; u < ((GRIDEXTENSION * 2) + 1); u++)
			{
				key[z] = std::make_tuple(GRID_x + gridoffsetx[u], GRID_y + gridoffsety[v], GRID_z + gridoffsetz[w]);

				//put empty texture on every quad
				//prevent flickering on textures
				mainQuad[z].set_texture(emptyTexture);

				z++;
			}
		}
	}
	
	GridFrustrum cache = gridFrustrum;
	gridFrustrum.clear();

	std::vector<CopyTuple> copyarray;

	std::vector<RefreshTuple> refresharray;

	z = 0;
	for (int w = 0; w < ((GRIDEXTENSION * 2) + 1); w++)
	{
		for (int v = 0; v < ((GRIDEXTENSION * 2) + 1); v++)
		{
			for (int u = 0; u < ((GRIDEXTENSION * 2) + 1); u++)
			{
				if (cache.count(key[z]) == 1)
				{
					copyarray.push_back(std::make_tuple(z, key[z]));
				}
				else
				{
					refresharray.push_back(std::make_tuple(z, GRID_x + gridoffsetx[u], GRID_y + gridoffsety[v], GRID_z + gridoffsetz[w]));
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


	for (int i = 0; i < copyarray.size(); i++)
	{
		CopyTuple copy = copyarray.at(i);

		CopyAndRefreshTexture(copy, cache);
	}

	for (int i = 0; i < refresharray.size(); i++)
	{
		RefreshTuple refresh = refresharray.at(i);
		
		_alltasks[std::get<0>(refresh)]->index = std::get<0>(refresh);
		_alltasks[std::get<0>(refresh)]->gridx = std::get<1>(refresh);
		_alltasks[std::get<0>(refresh)]->gridy = std::get<2>(refresh);
		_alltasks[std::get<0>(refresh)]->gridz = std::get<3>(refresh);

		//CleanTexture(gridTextureArray[std::get<0>(refresh)]);
		//gridTextureArray[std::get<0>(refresh)]->set_clear_color((0.0, 0.0, 0.0, 0.0));
		//gridTextureArray[std::get<0>(refresh)]->clear_image();

		GenericAsyncTask* refreshtask = new GenericAsyncTask("Refresh Grid", &renderParallelTextureTask, _alltasks[std::get<0>(refresh)]);
		
		refreshtask->set_task_chain("renderchain");

		taskMgr->add(refreshtask);
	}

}

void CopyTexture(PT(Texture) origin, PT(Texture) destination)
{
	int textsize = TEXTURESIZE * TEXTURESIZE * TEXTURESIZE * 4;

	PTA_uchar image = destination->modify_ram_image();
	//CPTA_uchar imageorg = origin->get_ram_image();
	PTA_uchar imageorg = origin->modify_ram_image();

	memcpy(image.p(), imageorg.p(), textsize);

	//int z = 0;
	//int y = 127;
	//int x = 0;

	//for (int i = 0; i < textsize; i += 4)
	//{
	//	image[i] = imageorg[i];  // BLUE COLOR
	//	image[i + 1] = imageorg[i + 1]; // GREEN COLOR
	//	image[i + 2] = imageorg[i + 2]; // RED COLOR
	//	image[i + 3] = imageorg[i + 3]; // ALPHA

	//	x++;

	//	if (x == TEXTURESIZE)
	//	{
	//		x = 0;
	//		y--;

	//		if (y == -1)
	//		{
	//			y = 127;
	//			z++;
	//		}
	//	}
	//}

	destination->reload();
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

	if (DENOISE == 1)
	{
		PT(Shader) denoise = Shader::load(Shader::ShaderLanguage::SL_GLSL, "shaders/denoise.vert", "shaders/denoise.frag");

		nodePath.set_shader(denoise);
	}
}

void InitShader(int index, NodePath nodePath)
{
	PT(TextureStage) ts = new TextureStage("ts");
	ts->set_mode(TextureStage::M_modulate);

	PT(Shader) myShader = Shader::load(Shader::ShaderLanguage::SL_GLSL, "shaders/shader.vert", "shaders/shader.frag");
	//nodePath.set_texture(ts, bunn);
	nodePath.set_shader_input("campos", camera.get_pos());
	nodePath.set_shader_input("params", LVector3f(GRID_SCALE, GRID_SCALE, INTERNALRES));
	nodePath.set_shader_input("target", mainWindow->get_render().get_relative_point(camera, LVector3f(0, 0, 1)));
	nodePath.set_shader(myShader);
	nodePath.set_transparency(TransparencyAttrib::Mode::M_alpha);

	mainQuad[index] = nodePath;
	nodePath.set_bin("fixed", (int) (index / ((GRIDEXTENSION * 2) + 1))); //I put 100 in the center quad WARN
}

void GenerateBillboard(int width, int height, WindowFramework * window, int index, bool useBuffer, NodePath parentNode, int centerx, int centery)
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
		InitShader(index, nodePath);

		nodePath.set_pos(centerx, 0, centery);
	}
	else
	{
		NodePath nodePath = window->get_pixel_2d().attach_new_node(node);
		InitShader(index, nodePath);
	}

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

	mybuffer = window->get_graphics_output()->make_texture_buffer("My Buffer", texsize, texsize);
	mytexture = mybuffer->get_texture();
	mybuffer->set_sort(-100);
	mycamera = new Camera("my camera");
	mycameraNP = window->get_render().attach_new_node(mycamera);
	region = mybuffer->make_display_region();
	region->set_camera(mycameraNP);
	myscene = NodePath("My Scene");
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
	
	GenerateMainBillboard(width, height, window, mytexture);

	for (int z = 0; z < ((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1)*((GRIDEXTENSION * 2) + 1); z++)
	{
		GenerateBillboard(width, height, window, z, true, myscene, -width / 2, height / 2);
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

	WindowProperties wp;
	framework.get_default_window_props(wp);
	int width = wp.get_x_size();
	int height = wp.get_y_size();

	window->enable_keyboard();

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

	// Get the camera and store it in a variable.
	camera = window->get_camera_group();
	camera.set_pos(0, 0, 1);
	//window->setup_trackball(); //move camera with mouse, errors while trying to move the camera from code directly

	initOffsetVectors();

	bunn = Render3dTexture(0, 0, 0);//dummy bunny for testing
	emptyTexture = Render3dTexture(0, 0, 0); 
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

	//RenderShadows(0, 0, 0);

	
	SceneGraphAnalyzer sga;
	sga.add_node(window->get_pixel_2d().node());
	sga.write(std::cerr);

	// Add our task.
	// If we specify custom data instead of NULL, it will be passed as the second argument
	// to the task function.
	AsyncTaskChain *chain = taskMgr->make_task_chain("changevdbgrid");
	chain->set_num_threads(1);
	chain->set_thread_priority(ThreadPriority::TP_urgent);

	renderChain = taskMgr->make_task_chain("renderchain");
	renderChain->set_num_threads(2);
	renderChain->set_thread_priority(ThreadPriority::TP_urgent);

	GenericAsyncTask* modifytask = new GenericAsyncTask("Modify Grid", &modifyGrid, nullptr);
	modifytask->set_task_chain("changevdbgrid");

	GenericAsyncTask* refreshtask = new GenericAsyncTask("Refresh Grid", &refreshGrid, nullptr);
	refreshtask->set_task_chain("changevdbgrid");

	
	taskMgr->add(new GenericAsyncTask("Camera Motion", &cameraMotionTask, nullptr));
	taskMgr->add(modifytask);
	taskMgr->add(refreshtask);

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
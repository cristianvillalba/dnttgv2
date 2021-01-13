// DnttgV2.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>

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

#include "VDBGrid.h"
#include "FilterManager.h"

#include "Dnntgv2.h"


#define WIDTH 200
#define HEIGHT 200
#define BUNNY 1
#define TEXTURESIZE 128

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

//Screen Projection Plane
PT(GeomVertexData) vdata;

//VDB class handler
VDBGrid*  grid;

//global 3d texture
PT(Texture) bunn;

//Main Quad
NodePath mainQuad;

//Window main
WindowFramework *mainWindow;

//Global camera pos;
float CAM_x = 0.0;
float CAM_y = 0.0;
float CAM_z = 1.0;

//Velocity vector/Acceleration
float VELX = 0.0;
float VELY = 0.0;
float ACCELERATION = 20.0;
float SPINVEL = 0;
float ANGLEDEGREES = 0;

//key states
bool FORWARD = false;
bool BACKWARDS = false;
bool LEFT = false;
bool RIGHT = false;
bool SPINL = false;
bool SPINR = false;



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
	grid->spawnSphere();

	refresh3dTexture();
}


PT(Texture) Render3dTexture()
{
	int texsize = TEXTURESIZE;

	PT(Texture)  bunn = new Texture("bunn");
	bunn->setup_3d_texture(texsize, texsize, texsize, Texture::ComponentType::T_float, Texture::Format::F_rgba8);

	for (int k = 0; k < texsize; k++) {
		PNMImage* pPNMImage = new PNMImage(texsize, texsize, 4);

		for (int i = 0; i < texsize; i++)
		{
			for (int j = 0; j < texsize; j++)
			{
				//This is good for index sample
				float x = (((float)i / texsize) - 0.5f) * 1000.0f;
				float y = (((float)j / texsize) - 0.5f) * 1000.0f + 250; //250 to put the bunny in the center of the render
				float z = (((float)k / texsize) - 0.5f) * 1000.0f;

				//this is good for world coordinates sample
				//float x = (((float)i / texsize) - 0.5f) * 100.0f;
				//float y = (((float)j / texsize) - 0.5f) * 100.0f + 25; //250 to put the bunny in the center of the render
				//float z = (((float)k / texsize) - 0.5f) * 100.0f;


				float r = 0.0f;
				float g = 0.0f;
				float b = 0.0f;

				//translate(500.0, 0, 0, &x, &y, &z);
				float data = grid->getValue(x, y, z);

				if (data > 0)
				{
					r = 1.0f ;
					g = 1.0f ;
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

	return bunn;
}

void refresh3dTexture()
{
	int textsize = TEXTURESIZE * TEXTURESIZE * TEXTURESIZE * 4;

	PTA_uchar image = bunn->modify_ram_image();

	int z = 0;
	int y = 127;
	int x = 0;

	for (int i = 0; i < textsize; i += 4)
	{
		float x0 = (((float)x / TEXTURESIZE) - 0.5f) * 1000.0f;
		float y0 = (((float)y / TEXTURESIZE) - 0.5f) * 1000.0f + 250; //250 to put the bunny in the center of the render
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
				y = 127;
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
	}

	if (BACKWARDS)
	{
		LVector3f fw = mainWindow->get_render().get_relative_point(camera, LVector3f(0, 0, 1.0));
		fw = fw - camera.get_pos();

		VELX = VELX + ACCELERATION * globalClock->get_dt() * fw.get_x();
		VELY = VELY + ACCELERATION * globalClock->get_dt() * fw.get_z();
	}

	if (LEFT)
	{
		LVector3f fw = mainWindow->get_render().get_relative_point(camera, LVector3f(-1.0, 0, 0.0));
		fw = fw - camera.get_pos();

		VELX = VELX + ACCELERATION * globalClock->get_dt() * fw.get_x();
		VELY = VELY + ACCELERATION * globalClock->get_dt() * fw.get_z();
	}

	if (RIGHT)
	{
		LVector3f fw = mainWindow->get_render().get_relative_point(camera, LVector3f(1.0, 0, 0.0));
		fw = fw - camera.get_pos();

		VELX = VELX + ACCELERATION * globalClock->get_dt() * fw.get_x();
		VELY = VELY + ACCELERATION * globalClock->get_dt() * fw.get_z();
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

	SPINVEL = SPINVEL * 0.9;

	if (abs(VELX) < 0.0001)
	{
		VELX = 0.0f;
	}

	if (abs(VELY) < 0.0001)
	{
		VELY = 0.0f;
	}

	if (abs(SPINVEL) < 0.0001)
	{
		SPINVEL = 0.0f;
	}

	CAM_x = CAM_x - VELX * globalClock->get_dt();
	CAM_z = CAM_z - VELY * globalClock->get_dt();
	ANGLEDEGREES = ANGLEDEGREES + SPINVEL * globalClock->get_dt();

	camera.set_pos(CAM_x, CAM_y,  CAM_z);
	camera.set_hpr(0, 0, ANGLEDEGREES);

	LVector3f lookAtDirection = mainWindow->get_render().get_relative_point(camera, LVector3f(0,0,1));

	mainQuad.set_shader_input("campos", camera.get_pos());
	mainQuad.set_shader_input("target", lookAtDirection);
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

	framework.define_key("space", "modify grid", spawnSphere, nullptr);

	// Get the camera and store it in a variable.
	camera = window->get_camera_group();
	camera.set_pos(0, 0, 1);
	//window->setup_trackball(); //move camera with mouse, errors while trying to move the camera from code directly

	//PT(Texture) bunn = Render3dTexture();
	bunn = Render3dTexture();

	//Procedurally generate a point in space
	vdata = new GeomVertexData("vertex", GeomVertexFormat::get_v3c4(), Geom::UH_static);

	vdata->set_num_rows(4);
	GeomVertexWriter vertex(vdata, "vertex");
	GeomVertexWriter color(vdata, "color");

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

	PT(Shader) myShader =  Shader::load(Shader::ShaderLanguage::SL_GLSL, "shaders/shader.vert", "shaders/shader.frag");
	nodePath.set_texture(bunn);
	nodePath.set_shader_input("campos", camera.get_pos());
	nodePath.set_shader_input("target", mainWindow->get_render().get_relative_point(camera, LVector3f(0, 0, 1)));
	nodePath.set_shader(myShader);

	mainQuad = nodePath;

	SceneGraphAnalyzer sga;
	sga.add_node(window->get_pixel_2d().node());
	sga.write(std::cerr);

	// Add our task.
	// If we specify custom data instead of NULL, it will be passed as the second argument
	// to the task function.
	//AsyncTaskChain *chain = taskMgr->make_task_chain("openvdbgrid");
	//chain->set_num_threads(4);
	//chain->set_thread_priority(ThreadPriority::TP_urgent);

	//GenericAsyncTask* rendertask = new GenericAsyncTask("Camera Task", &cameraTask, nullptr);
	//rendertask->set_task_chain("openvdbgrid");

	//taskMgr->add(rendertask);
	taskMgr->add(new GenericAsyncTask("Camera Motion", &cameraMotionTask, nullptr));

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
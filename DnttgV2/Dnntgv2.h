#pragma once
#include <tuple>
#include <queue>

using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::milli;
using std::random_device;
using std::sort;
using std::vector;

typedef std::tuple<int, int, int> KeyTriple;
typedef std::map<KeyTriple, unsigned char *> GridFrustrum;

typedef std::tuple<int, KeyTriple> CopyTuple;
typedef std::tuple<int, int, int, int> RefreshTuple;

typedef std::map<unsigned char *, bool> UsedTextures;

//structure to z depth ordering
struct Zorder {
	LVector3f pos;
	NodePath quad;
	NodePath quadNorm;
	int index;
};

//---Support structures
struct TaskArgs {
	int index;
	int debug;
	int gridx;
	int gridy;
	int gridz;
};

struct CopyArgs {
	int index;
	PT(Texture) data;
};

void MakeBunny(int argc, char *argv[]);
void MakeShadertoy(int argc, char *argv[]);
void print_results(const char *const tag,
	high_resolution_clock::time_point startTime,
	high_resolution_clock::time_point endTime);

void translate(float tx, float ty, float tz, float * txout, float * tyout, float * tzout);
void rotatex(float angle, float * txout, float * tyout, float * tzout);
void rotatey(float angle, float * txout, float * tyout, float * tzout);
void rotatez(float angle, float * txout, float * tyout, float * tzout);
void scale(float sf, float xf, float yf, float zf, float * txout, float * tyout, float * tzout);

void advanceCamera(const Event* eventPtr, void* dataPtr);
void advanceCameraUp(const Event* eventPtr, void* dataPtr);
void backCamera(const Event* eventPtr, void* dataPtr);
void backCameraUp(const Event* eventPtr, void* dataPtr);
void leftCamera(const Event* eventPtr, void* dataPtr);
void leftCameraUp(const Event* eventPtr, void* dataPtr);
void rightCamera(const Event* eventPtr, void* dataPtr);
void rightCameraUp(const Event* eventPtr, void* dataPtr);
void spinLCamera(const Event* eventPtr, void* dataPtr);
void spinLCameraUp(const Event* eventPtr, void* dataPtr);
void spinRCamera(const Event* eventPtr, void* dataPtr);
void spinRCameraUp(const Event* eventPtr, void* dataPtr);

void spawnSphere(const Event* eventPtr, void* dataPtr);

PT(Texture) Render3dTexture(int gridx, int griy, int gridz);
unsigned char * Render3dTextureAsArray(int gridx, int griy, int gridz);
PT(Texture) Render3dBigTexture();
void refresh3dTexture();
void refresh3dTexture(unsigned char * texture, int gridx, int gridy, int gridz);
void refresh3dTextureAsArray(unsigned char * texture, int gridx, int gridy, int gridz);
void callOpenGLSubImage(int posx, int posy, int posz, int debug);


//Headers
void initGridFrustrum();
void initOffsetVectors();
void initBigTexture();
void refreshGridFrustrum();
void GenerateMainBillboard(int w, int h, WindowFramework * window, PT(Texture) mtex);
void GeneratePrePassBillboard(int w, int h, WindowFramework * window, NodePath parent, int centerx, int centery);
void GenerateBillboard(int w, int h, WindowFramework * window, int index, bool useBuffer, NodePath parent, int centerx, int centery, int type);
void GenerateTextureBuffer(int w, int h, WindowFramework * window, NodePath s);
void InitShader(int index, NodePath node, int type);
void CopyTexture(PT(Texture) origin, PT(Texture) destination);
void CopyAndRefreshTexture(CopyTuple params, GridFrustrum cache, std::queue<TaskArgs> * copyQueue);
//void RefreshTexture(RefreshTuple params);
void RefreshTexture(KeyTriple params);
void CleanTexture(PT(Texture) origin);
unsigned char * CleanTextureArray();




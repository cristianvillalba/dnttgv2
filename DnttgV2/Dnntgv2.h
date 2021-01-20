#pragma once
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::milli;
using std::random_device;
using std::sort;
using std::vector;

typedef std::tuple<int, int, int> KeyTriple;
typedef std::map<KeyTriple, PT(Texture)> GridFrustrum;

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

PT(Texture) Render3dTexture();
void refresh3dTexture();

void initGridFrustrum();
void refreshGridFrustrum();
void GenerateMainBillboard(int w, int h, WindowFramework * window, PT(Texture) mtex);
void GenerateBillboard(int w, int h, WindowFramework * window, int index, bool useBuffer, NodePath parent, int centerx, int centery);
void GenerateTextureBuffer(int w, int h, WindowFramework * window, NodePath s);
void InitShader(int index, NodePath node);




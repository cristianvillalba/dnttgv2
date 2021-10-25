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
typedef std::map<KeyTriple, int> GridFrustrumPBO;

typedef std::tuple<int, KeyTriple> CopyTuple;
typedef std::tuple<int, int, int, int> RefreshTuple;

typedef std::map<unsigned char *, bool> UsedTextures;
typedef std::map<int, bool> UsedTexturesPBO;

//structure to z depth ordering
struct Zorder {
	LVector3f pos;
	NodePath quad;
	NodePath quadNorm;
	int index;
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
void refresh3dTexturePBO(int gridx, int gridy, int gridz, int pboindex, int saveptr);
int callOpenGLSubImage(int posx, int posy, int posz, int debug, int quad);
int callOpenGLSubImagePBO(int posx, int posy, int posz, int pboindex, int quad);



//Headers
void initGridFrustrum();
void initOffsetVectors();
void initBigTexture();
void initBigTexturePBO();
void refreshGridFrustrum();
void refreshGridFrustrumPBO();
void GenerateMainBillboard(int w, int h, WindowFramework * window, PT(Texture) mtex);
void GeneratePrePassBillboard(int w, int h, WindowFramework * window, NodePath parent, int centerx, int centery);
void GenerateBillboard(int w, int h, WindowFramework * window, int index, bool useBuffer, NodePath parent, int centerx, int centery, int type);
void GenerateTextureBuffer(int w, int h, WindowFramework * window, NodePath s);
void InitShader(int index, NodePath node, int type);
//void RefreshTexture(RefreshTuple params);
int RefreshTexture(KeyTriple params);
void CleanTexture(PT(Texture) origin);
unsigned char * CleanTextureArray();
void initCleanTexturePBO();
void exitCB();
int initPixelBufferObject();
void clearSharedMem();





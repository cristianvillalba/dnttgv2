#pragma once
//#include <thread.h>
#include <thread.h>
#include "VDBGrid.h"
#include "Dnntgv2.h"

class RefreshThread : public Thread
{
private:
	KeyTriple refresh;
	GridFrustrum gridFrustrum;
	std::queue<KeyTriple> * renderQueue;
	VDBGrid * grid;
	int texturesize;
	int boundingbox;

	void refresh3dTextureAsArray(unsigned char * texture, int gridx, int gridy, int gridz);
public:
	RefreshThread(KeyTriple refresh, GridFrustrum gridFrustrum, std::queue<KeyTriple> * renderQueue, VDBGrid * grid, int texturesize, int boundingbox);
	void thread_main();
};


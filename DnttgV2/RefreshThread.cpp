#include "RefreshThread.h"

RefreshThread::RefreshThread(KeyTriple ref, GridFrustrum gridF, std::queue<KeyTriple> * rQueue, 
	VDBGrid * gd, int texs, int bbox): Thread("refresh thread", "syncv1")
{
	refresh = ref;
	gridFrustrum = gridF;
	renderQueue = rQueue;
	grid = gd;
	texturesize = texs;
	boundingbox = bbox;
}

void RefreshThread::thread_main()
{
	std::cout << "refresh thread\n";

	int gridx = std::get<0>(refresh);
	int gridy = std::get<1>(refresh);
	int gridz = std::get<2>(refresh);


	if (gridFrustrum[refresh] != nullptr)
	{
		//Need to invert x and y axis
		refresh3dTextureAsArray(gridFrustrum[refresh], gridy, gridx, gridz);
		renderQueue->push(refresh);
		//return 0;
	}
	else
	{
		//return 1;
	}

	std::cout << "refresh thread - done!\n";
}

void RefreshThread::refresh3dTextureAsArray(unsigned char * texture, int gridx, int gridy, int gridz)
{
	int texsize = texturesize;

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

	//omp_set_dynamic(0);     // Explicitly disable dynamic teams
	//omp_set_num_threads(4); // Use 4 threads for all consecutive parallel regions
	//#pragma omp parallel for shared(tex, grid)
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

		//float data = grid->getValue(x - gridx * 1000, y - gridy * 1000, z - gridz * 1000);
		float data = grid->getValue(y - gridy * 1000, x - gridx * 1000, z - gridz * 1000);

		//if (data > 0)
		//if (data != BACKGROUNDVALUE)
		if (data != 1.5)
		{
			r = 255;
			g = 255;
		}

		if (boundingbox == 1)
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

		//#pragma omp critical
		tex[n] = r;
		tex[n + 1] = g;
		tex[n + 2] = b;
		tex[n + 3] = 255;
	}
}

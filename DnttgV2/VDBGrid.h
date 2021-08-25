#pragma once
#define NOMINMAX //to finally use movepoints -- conflict with MINMAX Macros
#define BACKGROUNDVALUE 666.0
#include <iostream>
#define _USE_MATH_DEFINES // for C++
#include <math.h>
#include <sstream>
#include <openvdb/points/PointCount.h>
#include <openvdb/points/PointMove.h>
#include <openvdb/openvdb.h>
#include <openvdb/io/Stream.h>
#include <openvdb/tools/ChangeBackground.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/ValueTransformer.h>


#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/Composite.h>

#include "pandaFramework.h"
#include "pandaSystem.h"

class VDBGrid
{
private:
	openvdb::tools::GridSampler<openvdb::FloatGrid::ConstAccessor, openvdb::tools::PointSampler>* fastSampler;
	//openvdb::tools::GridSampler<openvdb::points::PointDataGrid::ConstAccessor, openvdb::tools::PointSampler>* fastSamplerPoints;
	openvdb::FloatGrid::Ptr grid;
	//openvdb::points::PointDataGrid::Ptr points;
public:
	int initGrid();
	int initFastAccessor();
	
	float getValue(float x, float y, float z);
	//void movePoints();

	void spawnSphere(LVector3f pos, float radius);
	void spawnBox(LVector3f pos);
	
	template<class GridType>
	void makeBox(GridType& grid, const openvdb::Vec3f& c);
};


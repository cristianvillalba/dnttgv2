#include "VDBGrid.h"

using openvdb::GridPtrVecPtr;

// Define a local function that doubles the value to which the given
// value iterator points.
struct Local {
	static inline void op(const openvdb::FloatGrid::ValueAllIter& iter) {
		iter.setValue(*iter * 2);
	}
};

// This custom deformer is also used in the TestPointMove unit tests.
struct OffsetDeformer
{
	OffsetDeformer(const openvdb::Vec3d& _offset)
		: offset(_offset) { }
	template <typename LeafIterT>
	void reset(const LeafIterT&, size_t /*idx*/) { }
	template <typename IndexIterT>
	void apply(openvdb::Vec3d& position, const IndexIterT&) const
	{
		position += offset;
		//float x = rand();
		//float y = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		//float z = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		
		//openvdb::Vec3d randomvector(0, 0.1, x);
		//position += randomvector;

	}
	openvdb::Vec3d offset;
};

int VDBGrid::initGrid()
{
	std::cout << "Grid loading<-------------------------------------\n";

	// Initialize the OpenVDB library.  This must be called at least
// once per program and may safely be called multiple times.
	openvdb::initialize();

	/*
	// Create an empty floating-point grid with background value 0.
	openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create();
	std::cout << "Testing random access:" << std::endl;
	// Get an accessor for coordinate-based access to voxels.
	openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
	// Define a coordinate with large signed indices.
	openvdb::Coord xyz(1000, -200000000, 30000000);
	// Set the voxel value at (1000, -200000000, 30000000) to 1.
	accessor.setValue(xyz, 1.0);
	// Verify that the voxel value at (1000, -200000000, 30000000) is 1.
	std::cout << "Grid" << xyz << " = " << accessor.getValue(xyz) << std::endl;
	// Reset the coordinates to those of a different voxel.
	xyz.reset(1000, 200000000, -30000000);
	// Verify that the voxel value at (1000, 200000000, -30000000) is
	// the background value, 0.
	std::cout << "Grid" << xyz << " = " << accessor.getValue(xyz) << std::endl;
	// Set the voxel value at (1000, 200000000, -30000000) to 2.
	accessor.setValue(xyz, 2.0);
	// Set the voxels at the two extremes of the available coordinate space.
	// For 32-bit signed coordinates these are (-2147483648, -2147483648, -2147483648)
	// and (2147483647, 2147483647, 2147483647).
	accessor.setValue(openvdb::Coord::min(), 3.0f);
	accessor.setValue(openvdb::Coord::max(), 4.0f);
	std::cout << "Testing sequential access:" << std::endl;
	// Print all active ("on") voxels by means of an iterator.
	for (openvdb::FloatGrid::ValueOnCIter iter = grid->cbeginValueOn(); iter; ++iter) {
		std::cout << "Grid" << iter.getCoord() << " = " << *iter << std::endl;
	}
	*/

	// Create a VDB file object.
	openvdb::io::File file("bunny.vdb");
	//openvdb::io::File file("bunny_cloud.vdb");
	// Open the file.  This reads the file header, but not any grids.
	file.open();
	// Loop over all grids in the file and retrieve a shared pointer
	// to the one named "LevelSetSphere".  (This can also be done
	// more simply by calling file.readGrid("density").)
	openvdb::GridBase::Ptr baseGrid;
	for (openvdb::io::File::NameIterator nameIter = file.beginName();
		nameIter != file.endName(); ++nameIter)
	{
		// Read in only the grid we are interested in.
		//if (nameIter.gridName() == "density") {
		if (nameIter.gridName() == "ls_bunny") {
			baseGrid = file.readGrid(nameIter.gridName());
		}
		else {
			std::cout << "skipping grid " << nameIter.gridName() << std::endl;
		}
	}
	file.close();
	// From the example above, "LevelSetSphere" is known to be a FloatGrid,
	// so cast the generic grid pointer to a FloatGrid pointer.
	grid = openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid);
	
	//Empty grid
	//grid = openvdb::FloatGrid::create();
	//grid->setGridClass(openvdb::GRID_LEVEL_SET);
	
	openvdb::tools::changeBackground(grid->tree(), BACKGROUNDVALUE);

	//for debug purposes
	spawnSphere(LVector3f(0, 1000, 1000), 200.0, 1.0); //axis are inverted  looks like 1000 is middle, 
										      //and extension of grid is 1000 units
										      // also y axis is pointing negative on UP vector

	spawnSphere(LVector3f( 1000, 1000, 1000), 100.0, 1.0);
	spawnSphere(LVector3f( 2000, 1000, 1000), 50.0, 1.0);

	spawnBox(LVector3f(0, -500, 0));

	return 0;
}

int VDBGrid::initFastAccessor()
{
	// Request a value accessor for accelerated access.
	// (Because value accessors employ a cache, it is important to declare
	// one accessor per thread.)
	openvdb::FloatGrid::ConstAccessor accessor = grid->getConstAccessor();
	// Instantiate the GridSampler template on the accessor type and on
	// a box sampler for accelerated trilinear interpolation.
	
	fastSampler = new openvdb::tools::GridSampler<openvdb::FloatGrid::ConstAccessor, openvdb::tools::PointSampler>(accessor, grid->transform());
	
	// Request a value accessor for accelerated access.
	// (Because value accessors employ a cache, it is important to declare
	// one accessor per thread.)
	//openvdb::points::PointDataGrid::ConstAccessor accessor = points->getConstAccessor();
	// Instantiate the GridSampler template on the accessor type and on
	// a box sampler for accelerated trilinear interpolation.

	//fastSamplerPoints = new openvdb::tools::GridSampler<openvdb::points::PointDataGrid::ConstAccessor, openvdb::tools::PointSampler>(accessor, points->transform());

	return 0;
}

float VDBGrid::getValue(float x, float y, float z)
{
	// Compute the value of the grid at fractional coordinates in index space.
	openvdb::FloatGrid::ValueType indexValue = fastSampler->isSample(openvdb::Vec3R(x, y, z));

	return indexValue;

	//openvdb::points::PointDataGrid::ValueType indexValue = fastSamplerPoints->wsSample(openvdb::Vec3R(x, y, z));

	//return (float)indexValue;
}

void VDBGrid::spawnBox(LVector3f pos)
{
	//axis are inverted  looks like 1000 is middle
	openvdb::tools::changeBackground(grid->tree(), 1.5);

	openvdb::FloatGrid::Ptr boxGrid = openvdb::FloatGrid::create(/*background value=*/10.0);

	makeBox(*boxGrid, /*center=*/openvdb::Vec3f(pos.get_x(), pos.get_y(), pos.get_z()));

	for (openvdb::FloatGrid::ValueOffIter iter = boxGrid->beginValueOff(); iter; ++iter) {
		if (iter.getValue() < 0.0) {
			iter.setValue(1.0);
			iter.setValueOff();
		}
	}

	openvdb::tools::csgUnion(*grid, *boxGrid);

	//openvdb::tools::changeBackground(grid->tree(), 0.0);
}

//axis are inverted  looks like 1000 is middle, 
//and extension of grid is 1000 units
// also y axis is pointing negative on UP vector
void VDBGrid::spawnSphere(LVector3f pos, float radius, float voxelsize) //Axis are inverted!
{
	//axis are inverted  looks like 1000 is middle

	// Generate a level set grid.
	openvdb::FloatGrid::Ptr sphereGrid =
		openvdb::tools::createLevelSetSphere<openvdb::FloatGrid>(/*radius=*/radius,
			/*center=*/openvdb::Vec3f(pos.get_x()*voxelsize, pos.get_y()*voxelsize, pos.get_z()*voxelsize), /*voxel size=*/1.0 * voxelsize, 3.0, true);



	sphereGrid->setGridClass(openvdb::GRID_LEVEL_SET);

	// Convert the level set sphere to a narrow-band fog volume, in which
	// interior voxels have value 1, exterior voxels have value 0, and
	// narrow-band voxels have values varying linearly from 0 to 1.
	const float outside = sphereGrid->background();
	const float width = 2.0 * outside;
	// Visit and update all of the grid's active values, which correspond to
	// voxels on the narrow band.
	for (openvdb::FloatGrid::ValueOnIter iter = sphereGrid->beginValueOn(); iter; ++iter) {
		float dist = iter.getValue();
		iter.setValue((outside - dist) / width);
	}


	for (openvdb::FloatGrid::ValueOffIter iter = sphereGrid->beginValueOff(); iter; ++iter) {
		if (iter.getValue() < 0.0) {
			iter.setValue(1.0);
			iter.setValueOff();
		}
	}

	openvdb::tools::changeBackground(sphereGrid->tree(), BACKGROUNDVALUE);

	openvdb::tools::csgUnion(*grid, *sphereGrid);
}

template<class GridType>
void VDBGrid::makeBox(GridType& grid, const openvdb::Vec3f& c)
{
	using ValueT = typename GridType::ValueType;
	// Distance value for the constant region exterior to the narrow band
	const ValueT outside = grid.background();
	// Distance value for the constant region interior to the narrow band
	// (by convention, the signed distance is negative in the interior of
	// a level set)
	const ValueT inside = -outside;
	// Use the background value as the width in voxels of the narrow band.
	// (The narrow band is centered on the surface of the sphere, which
	// has distance 0.)
	int padding = int(openvdb::math::RoundUp(openvdb::math::Abs(outside)));
	// The bounding box of the narrow band is 2*dim voxels on a side.
	int dim = int(400.0 + padding);
	// Get a voxel accessor.
	typename GridType::Accessor accessor = grid.getAccessor();
	// Compute the signed distance from the surface of the sphere of each
	// voxel within the bounding box and insert the value into the grid
	// if it is smaller in magnitude than the background value.
	openvdb::Coord ijk;
	int &i = ijk[0], &j = ijk[1], &k = ijk[2];
	for (i = c[0] - dim; i < c[0] + dim; ++i) {
		const float x2 = openvdb::math::Pow2(i - c[0]);
		for (j = c[1] - dim; j < c[1] + dim; ++j) {
			const float x2y2 = openvdb::math::Pow2(j - c[1]) + x2;
			for (k = c[2] - dim; k < c[2] + dim; ++k) {
				// The distance from the sphere surface in voxels
				const float dist = 0.5;
				// Convert the floating-point distance to the grid's value type.
				ValueT val = ValueT(dist);
				// Only insert distances that are smaller in magnitude than
				// the background value.
				if (val < inside || outside < val) continue;
				// Set the distance for voxel (i,j,k).
				accessor.setValue(ijk, val);
			}
		}
	}
	// Propagate the outside/inside sign information from the narrow band
	// throughout the grid.
	openvdb::tools::signedFloodFill(grid.tree());
}

/*
void VDBGrid::movePoints()
{
	// Create an OffsetDeformer that moves the points downwards in Y by 10 world-space units.
	openvdb::Vec3d offset(0, -25, 0);
	OffsetDeformer deformer(offset);
	// Move the points using this deformer
	openvdb::points::movePoints(*points, deformer);
}
*/


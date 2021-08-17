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

	// Convert the level set sphere to a narrow-band fog volume, in which
	// interior voxels have value 1, exterior voxels have value 0, and
	// narrow-band voxels have values varying linearly from 0 to 1.
	
	//onst float outside = grid->background();
	//const float width = 2.0 * outside;
	
	// Visit and update all of the grid's active values, which correspond to
	// voxels on the narrow band.
	//for (openvdb::FloatGrid::ValueOnIter iter = grid->beginValueOn(); iter; ++iter) {
		//float dist = iter.getValue();
		//iter.setValue((outside - dist) / width);
	//}
	// Visit all of the grid's inactive tile and voxel values and update the values
	// that correspond to the interior region.
	for (openvdb::FloatGrid::ValueOffIter iter = grid->beginValueOff(); iter; ++iter) {
		if (iter.getValue() < 0.0) {
			iter.setValue(1.0);
			iter.setValueOff();
		}
	}
	// Set exterior voxels to 0.
	openvdb::tools::changeBackground(grid->tree(), 0.0);

	//check values from grid in file
	//for (openvdb::FloatGrid::ValueOnCIter iter = grid->cbeginValueOn(); iter; ++iter) {
	//	std::cout << "Grid" << iter.getCoord() << " = " << *iter << " value = " << iter.getValue() << std::endl;
	//}

	/*
	// Retrieve the number of leaf nodes in the grid.
	openvdb::Index leafCount = grid->tree().leafCount();
	// Use the topology to create a PointDataTree
	openvdb::points::PointDataTree::Ptr pointTree(
		new openvdb::points::PointDataTree(grid->tree(), 0, openvdb::TopologyCopy()));
	// Ensure all tiles have been voxelized
	pointTree->voxelizeActiveTiles();
	// Define the position type and codec using fixed-point 16-bit compression.
	using PositionAttribute = openvdb::points::TypedAttributeArray<openvdb::Vec3f,
		openvdb::points::FixedPointCodec<false>>;
	openvdb::NamePair positionType = PositionAttribute::attributeType();
	// Create a new Attribute Descriptor with position only
	openvdb::points::AttributeSet::Descriptor::Ptr descriptor(
		openvdb::points::AttributeSet::Descriptor::create(positionType));
	// Determine the number of points / voxel and points / leaf.
	openvdb::Index pointsPerVoxel = 1;
	openvdb::Index voxelsPerLeaf = openvdb::points::PointDataGrid::TreeType::LeafNodeType::SIZE;
	openvdb::Index pointsPerLeaf = pointsPerVoxel * voxelsPerLeaf;
	// Iterate over the leaf nodes in the point tree.
	for (auto leafIter = pointTree->beginLeaf(); leafIter; ++leafIter) {
		// Initialize the attributes using the descriptor and point count.
		leafIter->initializeAttributes(descriptor, pointsPerLeaf);
		// Initialize the voxel offsets
		openvdb::Index offset(0);
		for (openvdb::Index index = 0; index < voxelsPerLeaf; ++index) {
			offset += pointsPerVoxel;
			leafIter->setOffsetOn(index, offset);
		}
	}
	// Create the points grid.
	points =
		openvdb::points::PointDataGrid::create(pointTree);
	// Set the name of the grid.
	points->setName("Points");
	// Copy the transform from the sphere grid.
	points->setTransform(grid->transform().copy());
	// Randomize the point positions.
	std::mt19937 generator(0);
	std::uniform_real_distribution<> distribution(-0.5, 0.5);
	// Iterate over the leaf nodes in the point tree.
	for (auto leafIter = points->tree().beginLeaf(); leafIter; ++leafIter) {
		// Create an AttributeWriteHandle for position.
		// Note that the handle only requires the value type, not the codec.
		openvdb::points::AttributeArray& array = leafIter->attributeArray("P");
		openvdb::points::AttributeWriteHandle<openvdb::Vec3f> handle(array);
		// Iterate over the point indices in the leaf.
		for (auto indexIter = leafIter->beginIndexOn(); indexIter; ++indexIter) {
			// Compute a new random position (in the range -0.5 => 0.5).
			openvdb::Vec3f positionVoxelSpace(distribution(generator));
			// Set the position of this point.
			// As point positions are stored relative to the voxel center, it is
			// not necessary to convert these voxel space values into
			// world-space during this process.
			handle.set(*indexIter, positionVoxelSpace);
		}
	}
	// Verify the point count.
	openvdb::Index count = openvdb::points::pointCount(points->tree());
	std::cout << "LeafCount=" << leafCount << std::endl;
	std::cout << "PointCount=" << count << std::endl;
	*/


	//for debug purposes
	spawnSphere(LVector3f(0, 1000, 1000), 200.0); //axis are inverted  looks like 1000 is middle, 
										      //and extension of grid is 1000 units
										      // also y axis is pointing negative on UP vector

	spawnSphere(LVector3f( 1000, 1000, 1000), 100.0);
	spawnSphere(LVector3f( 2000, 1000, 1000), 50.0);

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

	openvdb::tools::changeBackground(grid->tree(), 0.0);
}

//axis are inverted  looks like 1000 is middle, 
//and extension of grid is 1000 units
// also y axis is pointing negative on UP vector
void VDBGrid::spawnSphere(LVector3f pos, float radius) //Axis are inverted!
{
	//axis are inverted  looks like 1000 is middle
	openvdb::tools::changeBackground(grid->tree(), 1.5);

	// Generate a level set grid.
	openvdb::FloatGrid::Ptr sphereGrid =
		openvdb::tools::createLevelSetSphere<openvdb::FloatGrid>(/*radius=*/radius,
			/*center=*/openvdb::Vec3f(pos.get_x(), pos.get_y(), pos.get_z()), /*voxel size=*/1.0);

	for (openvdb::FloatGrid::ValueOffIter iter = sphereGrid->beginValueOff(); iter; ++iter) {
		if (iter.getValue() < 0.0) {
			iter.setValue(1.0);
			iter.setValueOff();
		}
	}


	openvdb::tools::csgUnion(*grid, *sphereGrid);

	openvdb::tools::changeBackground(grid->tree(), 0.0);

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


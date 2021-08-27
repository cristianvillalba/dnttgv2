#version 430

// Set the number of invocations in the work group.
// In this case, we operate on the image in 16x16 pixel tiles.
layout (local_size_x = 16, local_size_y = 16) in;

// Declare the texture inputs
uniform readonly image3D fromTex;
uniform writeonly image3D toTex;

void main() {
  // Acquire the coordinates to the texel we are to process.
  ivec3 texelCoords = ivec3(gl_GlobalInvocationID.xyz);

  // Read the pixel from the first texture.
  vec4 pixel = imageLoad(fromTex, texelCoords);

  // Swap the red and green channels.
  //pixel.rg = pixel.gr;

  // Now write the modified pixel to the second texture.
  imageStore(toTex, texelCoords, pixel);
}
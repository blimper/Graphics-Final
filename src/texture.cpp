#include "texture.h"

#include <assert.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace CMU462 {

inline void uint8_to_float( float dst[4], unsigned char* src ) {
  uint8_t* src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
  dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8( unsigned char* dst, float src[4] ) {
  uint8_t* dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[0])));
  dst_uint8[1] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[1])));
  dst_uint8[2] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[2])));
  dst_uint8[3] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[3])));
}

void Sampler2DImp::generate_mips(Texture& tex, int startLevel) {

  // NOTE(sky): 
  // The starter code allocates the mip levels and generates a level 
  // map simply fills each level with a color that differs from its
  // neighbours'. The reference solution uses trilinear filtering
  // and it will only work when you have mipmaps.

  // Task 7: Implement this

  // check start level
  if ( startLevel >= tex.mipmap.size() ) {
    std::cerr << "Invalid start level"; 
  }

  // allocate sublevels
  int baseWidth  = tex.mipmap[startLevel].width;
  int baseHeight = tex.mipmap[startLevel].height;
  int numSubLevels = (int)(log2f( (float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  tex.mipmap.resize(startLevel + numSubLevels + 1);

  int width  = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel& level = tex.mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width  = max( 1, width  / 2); assert(width  > 0);
    height = max( 1, height / 2); assert(height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(4 * width * height);

  }
for(size_t level = 1; level < tex.mipmap.size(); ++level) {
    MipLevel& mip = tex.mipmap[level];
    for(int y = 0; y < mip.height; y++)
    {
      for(int x = 0; x < mip.width; x++)
      {
          float r_avg = 0;
          float g_avg = 0;
          float b_avg = 0;
          float a_avg = 0;  
          for(int i = 0; i < level; i++)
          {
              for(int j = 0; j < level; j++)
              {
                  r_avg += tex.mipmap[0].texels[4 * (x * level + i + (y * level + j) * mip.width * level)    ];
                  g_avg += tex.mipmap[0].texels[4 * (x * level + i + (y * level + j) * mip.width * level) + 1];
                  b_avg += tex.mipmap[0].texels[4 * (x * level + i + (y * level + j) * mip.width * level) + 2];
                  a_avg += tex.mipmap[0].texels[4 * (x * level + i + (y * level + j) * mip.width * level) + 3];
                }
            }

            r_avg = r_avg / (level * level);
            g_avg = g_avg / (level * level);
            b_avg = b_avg / (level * level);
            a_avg = a_avg / (level * level);

            mip.texels[4 * (x + y * mip.width)    ] = (uint8_t) (r_avg);
            mip.texels[4 * (x + y * mip.width) + 1] = (uint8_t) (g_avg);
            mip.texels[4 * (x + y * mip.width) + 2] = (uint8_t) (b_avg);
            mip.texels[4 * (x + y * mip.width) + 3] = (uint8_t) (a_avg);


          
      }
  }
}
/*
  // fill all 0 sub levels with interchanging colors
  Color colors[3] = { Color(1,0,0,1), Color(0,1,0,1), Color(0,0,1,1) };
  for(size_t i = 1; i < tex.mipmap.size(); ++i) {

    Color c = colors[i % 3];
    MipLevel& mip = tex.mipmap[i];

    for(size_t i = 0; i < 4 * mip.width * mip.height; i += 4) {
      float_to_uint8( &mip.texels[i], &c.r );
    }
  }
  */
}

Color Sampler2DImp::sample_nearest(Texture& tex, 
                                   float u, float v, 
                                   int level) {

  // Task ?: Implement nearest neighbour interpolation
     MipLevel *miplevel = &tex.mipmap[level];
     int offset = (int)(floor(miplevel->width * u)) + (int)(floor(miplevel->height * v) * miplevel->width);
     offset = offset * 4;
     return Color((float)miplevel->texels[offset]/255, 
             (float)miplevel->texels[offset+1]/255, 
             (float)miplevel->texels[offset+2]/255, 
             (float)miplevel->texels[offset+3]/255);
  // return magenta for invalid level
     //return Color(1,0,1,1);

}

float colorWeight(MipLevel *miplevel, int index, int offset, float diffx, float diffy)
{
  int index1 = index + 4;
  int index2 = index + miplevel->width * 4;
  int index3 = index2 + 4;
  float rdiffx = 1 - diffx;
  float rdiffy = 1 - diffy;

  float r = (float)miplevel->texels[index + offset] * diffx;
  r += (float)miplevel->texels[index1 + offset] * rdiffx;
  float r2 = (float)miplevel->texels[index2 + offset] * diffx;
  r2 += (float)miplevel->texels[index3 + offset] * rdiffx;
  r = r * diffy + r2 * rdiffy;
  return r;
}
Color Sampler2DImp::sample_bilinear(Texture& tex, 
                                    float u, float v, 
                                    int level) {
  
  // Task ?: Implement bilinear filtering
  MipLevel *miplevel = &tex.mipmap[level];
  int texx = (int)floor(miplevel->width * u - 0.5);
  int texy = (int)floor(miplevel->height * v - 0.5);
  float diffx = miplevel->width * u - ((float)texx + 0.5);
  float diffy = miplevel->height * v - ((float)texy + 0.5);
  float rdiffx = 1 - diffx;
  float rdiffy = 1 - diffy;
  
  int index = (int)(floor(miplevel->width * u)) + (int)(floor(miplevel->height * v) * miplevel->width);
  index = index * 4;



  float r = colorWeight(miplevel, index, 0, diffx, diffy);

  float g = colorWeight(miplevel, index, 1, diffx, diffy);

  float b = colorWeight(miplevel, index, 2, diffx, diffy);

  float a = colorWeight(miplevel, index, 3, diffx, diffy);

  return Color((float)r/255, (float)g/255, (float)b/255, (float)a/255);
  // return magenta for invalid level
  //return Color(1,0,1,1);

}

Color Sampler2DImp::sample_trilinear(Texture& tex, 
                                     float u, float v, 
                                     float u_scale, float v_scale) {

  // Task 8: Implement trilinear filtering
    int baseWidth  = tex.mipmap[0].width;
    int baseHeight = tex.mipmap[0].height;
    int numSubLevels = (int)(log2f( (float)max(baseWidth, baseHeight)));

    numSubLevels = min(numSubLevels, kMaxMipLevels - 1);

    float dudx = 1/u_scale * tex.width;
    float dvdy = 1/v_scale * tex.height;
    float L = max(dudx, dvdy);
    float d = log(L);
    if (d <= 0)
        return sample_bilinear(tex, u, v);
    else if (d >= numSubLevels)
        return sample_bilinear(tex, u, v, numSubLevels);
    else
    {
        int level = (int)floor(d);
        float distance = d - level;
        Color c1 = sample_bilinear(tex, u, v, level + 1);
        Color c2 = sample_bilinear(tex, u, v, level);
        return Color(c1.r * distance + c2.r * (1-distance), 
                c1.g * distance + c2.g * (1-distance), 
                c1.b * distance + c2.b * (1-distance), 
                c1.a * distance + c2.a * (1-distance));
    }

}

} // namespace CMU462

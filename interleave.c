#include <stdio.h>
#include <stdint.h>

int main(int argc, char **argv)
{
  if(argc != 4)
  {
    printf("Usage: %s infile1 infile2 outfile\n", *argv);
    printf("  interleaves two files");
  }
  else
  {
    FILE *f1 = fopen(argv[1], "rb");
    FILE *f2 = fopen(argv[2], "rb");
    FILE *of = fopen(argv[3], "wb");
    if(!f1 || !f2 || !of)
      err(1, "fopen");
    uint16_t i;
    while(!feof(f1) && !feof(f2))
    {
      fread(&i, 1, 2, f1);
      fwrite(&i, 1, 2, of);
      fread(&i, 1, 2, f2);
      fwrite(&i, 1, 2, of);
    }
    fclose(of);
  }
  return 0;
}

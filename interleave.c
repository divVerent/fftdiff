#include <stdio.h>
#include <stdint.h>

int main(int argc, char **argv)
{
  FILE *f1 = fopen(argv[1], "rb");
  FILE *f2 = fopen(argv[2], "rb");
  FILE *of = fopen(argv[3], "wb");
  uint16_t i;
  while(!feof(f1) && !feof(f2))
  {
    fread(&i, 1, 2, f1);
    fwrite(&i, 1, 2, of);
    fread(&i, 1, 2, f2);
    fwrite(&i, 1, 2, of);
  }
  fclose(of);
  return 0;
}

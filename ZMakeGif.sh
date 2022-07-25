if [ "$#" -ne 1 ]; then
  echo "Usage: $0 PngFilePrefix" >&2
  exit 1
fi

ffmpeg -y -f image2 -framerate 5 \
-i $1%003d.png \
-pix_fmt rgb8 \
-vf scale=1280:-1 \
$1.mp4

ffmpeg -y -i $1.mp4 -vf "scale=1280:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" $1.gif
rm $1.mp4
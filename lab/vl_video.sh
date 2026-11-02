#!/bin/bash
# VIDEO ATOMO: assemble the mp4s from the 1080p frames (1.5 s still at the start, 3 s at the end)
SCR=${VL_SCR:-/tmp/claude-0/-home-claude-particle-sim/a13a264e-04a0-55c0-8c14-0e293e6cf6fd/scratchpad/vl}
F="$SCR/frames_1920/%05d.png"
ffmpeg -y -loglevel error -framerate 24 -i "$F" -vf "tpad=start_duration=1.5:start_mode=clone:stop_duration=3:stop_mode=clone" \
  -c:v libx264 -preset slow -crf 18 -pix_fmt yuv420p -movflags +faststart "$SCR/video_atomo_1080p.mp4"
ffmpeg -y -loglevel error -framerate 24 -i "$F" -vf "tpad=start_duration=1.5:start_mode=clone:stop_duration=3:stop_mode=clone,scale=1280:720:flags=lanczos" \
  -c:v libx264 -preset slow -crf 20 -pix_fmt yuv420p -movflags +faststart "$SCR/video_atomo_720p.mp4"
ls -la "$SCR"/video_atomo_*.mp4

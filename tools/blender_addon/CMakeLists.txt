cmake_minimum_required(VERSION 3.0)

if(BUILD_BLENDER)
  find_package(Blender REQUIRED)

  # message(STATUS "I'm currently looking at the blender python installation")
  install(DIRECTORY . DESTINATION ${BLENDER_ADDON_PATH}/addons/gamer_addon FILES_MATCHING PATTERN "*.py")
endif(BUILD_BLENDER)
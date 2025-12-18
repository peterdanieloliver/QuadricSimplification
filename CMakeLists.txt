cmake_minimum_required(VERSION 3.5)
project(QuadricSimplification VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(Qt5_DIR "${Qt5_PATH}/lib/cmake/Qt5")
set(Qt5_BIN "${Qt5_PATH}/bin")
find_package(Qt5 REQUIRED COMPONENTS Core Widgets Gui)

set(SRC src)
set(SOURCES
    ${SRC}/main.cpp
    ${SRC}/dialogs.hpp
    ${SRC}/glwidget.hpp
    ${SRC}/elements.hpp
    ${SRC}/icMatrix.hpp
    ${SRC}/icVector.hpp
    ${SRC}/mainwindow.hpp
    ${SRC}/messagestream.hpp
    ${SRC}/ply.cpp
    ${SRC}/polyhedron.hpp
    ${SRC}/scene.hpp
    ${SRC}/trackball.hpp
    ${SRC}/window.hpp
)
add_executable(QuadricSimplification ${SOURCES})

target_link_libraries(QuadricSimplification PRIVATE Qt5::Core Qt5::Widgets Qt5::Gui)

target_include_directories(QuadricSimplification PRIVATE include stb)

# Copy Qt5 DLLs to the executable directory on Windows
if(WIN32)

    set(QT_DLLS
        Qt5Core
        Qt5Widgets
        Qt5Gui
    )

    foreach(QT_DLL ${QT_DLLS})
        add_custom_command(
            TARGET QuadricSimplification POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<$<CONFIG:Debug>:"${Qt5_BIN}/${QT_DLL}d.dll">
                $<$<NOT:$<CONFIG:Debug>>:"${Qt5_BIN}/${QT_DLL}.dll">
                $<TARGET_FILE_DIR:QuadricSimplification>
        )
    endforeach()

     # Copy Qt platform plugins (e.g., qwindows.dll)
    add_custom_command(
        TARGET QuadricSimplification POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory
            $<TARGET_FILE_DIR:QuadricSimplification>/platforms
    )

    add_custom_command(
        TARGET QuadricSimplification POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<$<CONFIG:Debug>:"${Qt5_PATH}/plugins/platforms/qwindowsd.dll">
            $<$<NOT:$<CONFIG:Debug>>:"${Qt5_PATH}/plugins/platforms/qwindows.dll">
            $<TARGET_FILE_DIR:QuadricSimplification>/platforms
    )

    # Copy the entire platforms directory
    # add_custom_command(
    #     TARGET QuadricSimplification POST_BUILD
    #     COMMAND ${CMAKE_COMMAND} -E copy_directory
    #         "${Qt5_PATH}/plugins/platforms"
    #         $<TARGET_FILE_DIR:QuadricSimplification>/platforms
    # )

endif()

***************************************
Using GAMer in Blender to Refine a Mesh
***************************************

Introduction
---------------------------------------------

In this tutorial, we will use GAMer to improve the quality of a surface mesh.

Our model of interest today is the murine cardiac myocyte calcium release unit.
You can get the segmentation from the Cell Centered Database (CCDB), a
collection of 2D, 3D, and 4D cellular and subcellular data derived from light
and electron microscopy. The data we will be using today is from entry 3603.
For the purposes of this tutorial, we have selected a nicely segmented section
for model building. IMOD was used to generate a preliminary surface mesh in obj
format.

Preliminary surface meshes from segmentations are typically of poor quality for
computation. A poor quality mesh is defined as one which will lead to large
errors in computation or inefficiency. For triangle surface meshes, the quality
of the mesh can be qualitatively judged by evenness of the edge and node
distribution. Quantitatively, mesh quality is measured by mesh topology and
proportion of high aspect triangles (triangles with extremely large and small
angles).

Import OBJ File into Blender
---------------------------------------------

- Download the `obj file`_.

.. _obj file: https://github.ctlee/gamer/examples/tt-sr-mit.obj

- Now we can open up this ".obj" file in blender.

  - Start Blender by typing ``blender`` at the command line.

  - Change the Rotation Mode by going to **User Preferences > Input tab >
    track-ball**

  - Import the data file with **File > Import > Wavefront (.obj)** and select
    **tt-sr-mit.obj**

  .. image:: ../images/import_obj.png


Preliminary Work on Imported Mesh
---------------------------------------------

- The object is not centered around the origin. To bring it into view, do the
  following.
    - Make sure you have **all four objects selected** by holding **shift-LMB** and **clicking**
      (i.e. left-mouse-button) each object, selected objects will have a
      orange circle around their respective triangle.
    - Have the cursor be inside of the **3D view window** and press the **Period** on
      the **Numpad**.


  .. image:: ../images/3dview_object_selected.png

- You may notice that parts of the model are getting truncated by the clipping
  plane. To remove the visual artifacts, we perform the following:

  - Increase the distance of the far clipping plane

    - Open the **Properties** panel by  having the cursor in the **3D view window** and
      then hit  **n**
    - Navigate to the **View** subpanel
    - Under **Clip**, change **End** to **2000**.

    .. image:: ../images/end_2000.png

  - For geometric models, it’s often useful to change to the orthographic view
    [**numpad 5**].

- We need the model to be in one volumetric domain so let’s join it. In the
  **Outliner** hold **shift-LMB** and **click** (i.e. left-mouse-button)  each object
  with **obj1_T-Tub_1** selected last which will appear white.

  .. image:: ../images/objects_selected.png

- Then to join them  have the cursor in the **3D view window** and press (**Ctrl+j**), if
  done correctly then the four objects should combine into **obj1_T-Tub_1**.

  .. image:: ../images/combined_object_selected.png

- The mesh is currently rendered as a solid material. While this is great for
  the purposes of animation and visualization, we care about the distribution
  of nodes and edges. To show the edges in the Viewport, go to the **Object
  Properties** tab and select **Wire** and **Draw all Edges** in the
  **Display** subpanel.

  .. image:: ../images/wire_draw_all_edges.png

- To simplify future manipulation let’s center the model about the origin.

  - select **Object** near the bottom left of the 3D window, then select **transform**,
    then select **Geometry to Origin**.

  .. image:: ../images/object_origin.png

  - Set the origin again just as before with the **Period** on the **Numpad** then to set
    the focus to the front of the object press **1** on the **Numpad**. Your **3D View window**
    should now look something like the following.

    .. image:: ../images/front_face_object.png

- Let’s now align the model so that the long axis is horizontal.

  - Rotate about the y-axis by 45 degrees to line up the model horizontally, by
    hitting **r**, **y**, and **45**.

    .. image:: ../images/front_face_object_horiz.png

  - Save this state as the object’s default rotation and scale via one of two ways.

    - Select **Object** near the bottom left of the window, then select **Apply**,
      then select **Rotation and Scale**.

      .. image:: ../images/object_apply_rotationscale.png

    - Or you can press **Ctrl+a** and then select **Rotation and Scale**.

      .. image:: ../images/ctrl_rotationscale.png

- CHECKPOINT: Let’s save our work now as: **tt-sr-mit.imp_obj.blend**. Note
  that if something goes awry, you can always close Blender and reopen at this
  checkpoint!

Analyze Mesh, Clean-up, and Repeat
---------------------------------------------

- We can use the mesh analyzer function of CellBlender to inspect the mesh for
  suitability for computational analysis. I.e., that there are no unexpected
  holes or bizarre topologies.

  - Click on the **CellBlender-Mesh Analysis** drop down located just above the
    **CellBlender** dropdown, then press the **Analyze Mesh** button.You should
    see that it is **not watertight and non-manifold**. Now we know that there is
    a hole in the mesh somewhere, rendering it non-watertight. Similarly there are some
    topological issues indicated by non-manifold topology. Manifold geometry is
    essentially geometry which can exist in the real world. For some pragmatic
    examples of non-manifold geometry please consult the following
    stackexchange.

    .. image:: ../images/analyze_mesh.png

- Let’s start by cleaning up regions of non-manifold topology.

  - First engage **Edit Mode** [**Tab**] and while having the cursor in the **3D view window**
    deselect everything by pressing [**a**].
  - Hit **Ctrl-Tab** and select **Vertex** select mode.

    .. image:: ../images/vertex_select.png

  - Click **Select** near the bottom left of the window, then click **Select All By Trait**,
    then click **Non Manifold**.

    .. image:: ../images/select_selectbytrait_nonmanifold.png
  - Or you could press [**Shift+Ctrl+Alt+m**] as a shortcut.

  - This highlights all the regions of **non-manifold topologies**.

    .. image:: ../images/non_manifold.png

- Conveniently non-manifoldness is a problem in the animation industry (it
  tends to cause problems with raytracing among other things). Thus, Blender
  has some built-in tools to help resolve non-manifoldness.

  - First, Select All by pressing [**a**] with the cursor in the **3D view window**, then near the botttom left of the 3D window
    select **Mesh**, then **Clean up**, then **Degenerate** and finally **Dissolve**. This function will take care
    of several cases of bad geometry: edges with no length, faces with no area, or face corners with no area. It does
    so by deleting vertices and edges it thinks don’t make sense.

    .. image:: ../images/degenerate_dissolve.png

  - This will leave some holes in the mesh. We can automatically fill the holes
    by again selecting **Mesh** near the botttom left of the 3D window, then **Clean up**, then **Fill Holes**.

    .. image:: ../images/fill_holes.png

  - Let’s now check how many issues we have resolved. Deselect everything by pressing [**a**] with the cursor in the
    **3D window** again and then near the botttom left of the 3D window  click **Select**, then **Select All By Trait**,
    then **Non Manifold**. Or we could use [**Shift+Ctrl+Alt+m**] as a shortcut.

  - We see that the mesh has been substantially improved but is not perfect yet.

    .. image:: ../images/almost_manifold.png

- We can zoom in on the selected region by again having the cursor in the 3D window and then on the **Numpad** select the
  **Period**.

  - Let’s delete the dangling vertex. First Deselect everything [**a**] then
    select the culprit vertex [**RMB click**] (**Note** this can be difficult to find so make sure you have the view **outside** the
    object and **not inside**) and delete [**x**] and choose Vertices.

    - Normal view of the culprit vertx

      .. image:: ../images/culprit_vertex.png

    - Close up of the culprit vertex

      .. image:: ../images/culprit_vertex_zoom.png

- Once again let’s take a look to see if there are any residual problems. In **Edit Mode** [**Tab**], click **Select**,
  then **Select All By Trait**, then **Non Manifold**. Or we could use [**Shift+Ctrl+Alt+m**] as a shortcut. At this
  point your mesh should have no more issues.
- Recall that the degenerate dissolve function deleted some vertices and edges.
  In some cases, when the holes are filled, the polygons may no longer be
  triangular.

  - To re-triangulate, select everything [**a**] and choose **Mesh**, then **Faces**, then **Triangulate**. Or [**Ctrl+t**]
    as a shortcut.

    .. image:: ../images/mesh_faces_triangle.png

- Our mesh is starting to look pretty good! Let’s re-run mesh analyzer

  - Return to **Object Mode** [**Tab**] or by pressing the list by the bottom of the 3d window.

    .. image:: ../images/tabbutton.png
    .. image:: ../images/tabbutton_objectmode.png

  - Rerun mesh analysis: click the drop down **CellBlender-Mesh Analysis**, then **Analyze Mesh**. We now
    have a **Watertight** and **Manifold** mesh but we have **Inward Facing normals**. This
    means that everything is good except the mesh is **inside out**!

    .. image:: ../images/analyze_mesh_fixed.png

- To reset the orientation of the faces, we need to recalculate the normals.

  - Return to **Edit Mode** [**Tab**].
  - Select **Mesh**, then **Normals**, then **Recalculate Outside** or you could use [**Ctrl+n**] as a shortcut.

    .. image:: ../images/mesh_normals_recalculate_outside.png

  - Return to to **Object Mode** [**Tab**], run mesh analyzer again. We now we have
    good geometry to start with. Be sure to note the **surface area** and **volume**.

    .. image:: ../images/analyze_mesh_area_volume.png

- CHECKPOINT: Save your progress to: **tt-sr-mit.clean.blend**.

Using GAMer
---------------------------------------------

- We are now ready to begin surface mesh refinement with GAMer.

  - Go to the **GAMer** tab on the left side of Blender.
  - Click on the **Surface Mesh Improvement** button to show this subpanel.

    .. image:: ../images/surface_mesh_improve.png

  - The subpanel provides several functions as follows:

    - **Coarse Dense Tris**: reduces the number of triangles in densely
      triangulated portions of the mesh.
    - **Coarse Flat Tris**: reduces the number of triangles in flat regions of
      the mesh.
    - **Smooth Tris**: improves the aspect ratio of triangles by maximizing
      angles. It does so by flipping edges moving vertices based on angle and
      the local structure tensor.
    - **Normal Smooth Surf**: smooths surface roughness using a
      feature-preserving normal averaging algorithm.

  - In **Object Mode** [**Tab**] with the model selected, perform the following
    operations in order. After each step the approximate number of vertices
    remaining is given.

    - **Smooth Tris**: Max_Min = 15, S_Iter = 10 (~73K vertices)

      .. image:: ../images/smooth_tris_changes.png

    - **Coarse Dense Tris**: CD_R, 1; CD_Iter, 5 (~37K vertices)

      .. image:: ../images/coarse_dense_tris_changes.png

    - **Smooth Tris**: Max_Min, 15; S_Iter, 10

      .. image:: ../images/smooth_tris_changes.png

    - **Coarse Dense Tris**: CD_R, 0.5; CD_Iter, 5 (~28K vertices)

      .. image:: ../images/coarse_dense_tris_decrement.png

    - **Smooth Tris**: Max_Min, 20; S_Iter, 20

      .. image:: ../images/smooth_tris_increment.png

    - Click **Normal Smooth Surf** twice

      .. image:: ../images/normal_smooth_surf_twice.png

  - While in **Object Mode** [**Tab**], click **CellBlender**, then **CellBlender-Mesh Analyzer**, then **Mesh Analyzer**.
    Note the slightly smaller **surface area** but similar **volume**.

      .. image:: ../images/analyze_mesh_area_volume_change.png

- CHECKPOINT: Save your progress to: **tt-sr-mit.gamer_proc_1.blend**

Add Boundary Box
---------------------------------------------

- Now that we have a reasonable surface mesh of our features, we want to place
  a boundary box around the features to represent the cytosol.

  - First we center the 3D cursor to the center. We will next add a cube at the
    position of the 3D cursor. In **Object Mode**, select **Object**, then **Snap**,
    then **Cursor to Center** or you could use [**Shift+s** and select **Cursor to Center**] as a shortcut.

    .. image:: ../images/object_snap_cursorcenter.png

  - We will next add a cube at the position of the 3D cursor. Add a cube mesh
    object, select **Add**, then **Mesh**, then **Cube**. Or you could use [**Shift+a** and select **Mesh**, then **Cube**]
    as a shortcut.

    .. image:: ../images/add_mesh_cube.png

  - Switch to **Edit mode** [**Tab**], let’s scale and translate the bounding box to where we want it. Recall that
    the **Properties** panel can be summoned with [**n**].

    - **Location** (-40, 15, 30)
    - **Scale** (275, 130, 220)

  .. image:: ../images/add_cube.png

- The cube is currently a quadrilateral mesh. We need to convert to a
  triangular mesh.

  - Switch to **Edit Mode** [**Tab**].
  - To capture detailed features we will need additional triangles. With the
    cube selected, select **Mesh**, then **Edges**, then **Subdivide** a total of six times. Or you could use [**w** and
    select **Subdivide**] as a shortcut.

  .. image:: ../images/mesh_edges_subdivide.png

  - Triangulate by selecting **Mesh**, then **Faces**, then **Triangulate Faces**. Or you could use [**Ctrl+t**] as a shortcut.
  - Return to **Object Mode** [**Tab**].

  .. image:: ../images/subdivide_cube.png

- CHECKPOINT: Save your progress to: **tt-sr-mit.with_cube.blend**

Using Boolean Modifier
---------------------------------------------

- To get the surface representation of the cytosolic volume, we must subtract
  our features from our cube mesh.

  - While in **Object Mode** [**Tab**], go to the **Modifier** tab of the
    **Properties Panel** and hit **Add Modifier**, **Generate: Boolean**,
    **Operation: Difference**, Object: **obj1_T-Tub_1** and **Apply** the
    modifier.
  - In the **Outliner** click on the eye to hide **obj1_T-tub_1**.
  - With the cube selected, apply the current rotation and scale transform.
    Select **Object**, then **Apply**, **Rotation and Scale**, or use [**Ctrl+a** and select
    **Rotation and Scale**]
  - Apply the current location transform. Select **Object**, then **Apply**, then **Location** or use
    [**Ctrl+a, Location**].
  - If you would like to show the edges, go to the **Object Properties** and
    select **Wire** and **Draw all Edges**.

  .. image:: ../images/add_boolean.png

- CHECKPOINT: Save your progress to: **tt-sr-mit.boolean.blend**


Refine Cube with GAMer
---------------------------------------------

- Once again, we have a surface mesh to refine.

  - First, in **Edit Mode** [**Tab**], switch to **Vertex** select mode.
  - Deselect everything [**a**].
  - Next, we can click **Select**, then **Select All By Trait**, then **Non Manifold**, or
    [**Shift+Ctrl+Alt+m**]. Nothing should be selected. If there are some
    issues, try performing **Degenerate Dissolve** followed by **Fill Holes**.
  - Return to **Object Mode** [**Tab**], and run **Mesh Analyzer**. We find
    that the mesh is not triangulated.

- We can triangulate as before:

  - In **Edit Mode** **Tab**, Select All [**a**] , then select **Mesh**, then **Faces**, then **Triangulate Faces** or [**Ctrl+t**]
  - Return to **Object Mode** [**Tab**], and run **Mesh Analyzer**. We have a good geometry to start refining.

- CHECKPOINT: Save your progress to: **tt-sr-mit.boolean_clean.blend**
- Let’s begin surface refinement using GAMer

  - In **Object Mode** [**Tab**] with the cube selected, perform the following
    operations in order. After each step the approximate number of vertices
    remaining is given.

    - **Smooth Tris**: Max_Min = 15, S_Iter = 10 (~70K vertices)
    - **Coarse Dense Tris**: CD_R = 0.75, CD_Iter = 10 (~57K vertices)
    - **Coarse Flat Tris**: CF_Rate = 0.016 (~44K vertices)
    - **Smooth Tris**: Max_Min = 15; S_Iter = 10
    - **Coarse Dense Tris**: CD_R = 0.1, CD_Iter = 10 (~42K vertices)
    - **Smooth Tris**: Max_Min = 20; S_Iter = 20
    - **Normal Smooth Surf** twice

  - In **Object Mode** [**Tab**], run **Mesh Analyzer**. Note the slightly
    smaller surface area but similar volume.

- CHECKPOINT: Save your progress to: **tt-sr-mit.gamer_proc_2.blend** Now we're
  ready to add boundaries and associated boundary markers to the mesh!


Adding Cytosolic Boundary
---------------------------------------------

- Return to the **GAMer** tab and choose the **Boundary Marking** tool

  - Add a new boundary (**+** button). By clicking on the color swatch, you can
    select the color you wish to represent the **Cytosol**. The color only
    serves as a visual aid to help you mark. Set the color to green.
  - Change the name of the boundary to **Cytosol**.

    .. image:: ../images/boundary_marking_cyto.png

  - Enter **Edit Mode** [**Tab**] and choose **Face** select mode and begin
    selecting all faces of the cytosol. Clicking each face is very arduous! For
    larger surfaces, you may elect to select using the **Circle Select** tool
    [**c**] or the **Border Select** tool [**b**]. Use "Assign" to assign
    selected faces to boundary. You can assign as you go or all together at the
    end. Note, it can sometimes be very helpful to hide all selected faces
    using [**h**], or hide all unselected faces using [**Shift+h**]. You can
    unhide everything using [**Alt+h**]. In the next steps, we'll be using the
    the **Border Select** tool [**b**].
  - Turn off the option: **Limit selection to visible**.
  - **Front-View** [**numpad 1**].
  - Select faces of **Cytosol**. Use **Border Select** tool [**b**] to select
    the profile of each side.
  - **Top-View** [**numpad 7**].
  - Select additional faces of **Cytosol**. Use **Border Select** tool [**b**]
    to select the profile of remaining sides.
  - Hide all unselected [**Shift+h**]. You may notice that some triangles from
    internal features may have been selected. We will fix this next by
    selecting linked triangles.
  - Deselect all [**a**]
  - Select one triangle, click [**RMB**].
  - Select Linked [**Ctrl+l**]
  - Hide All Deselected [**Shift+h**]
  - Use "Assign" to assign selected faces to boundary.
  - Turn on option: “Limit selection to visible”.
  - Unhide All [**Alt+h**]
  - Deselect all [**a**]

- CHECKPOINT: Save your progress to: **tt-sr-mit.cytosol.blend**


Adding Other Boundaries
---------------------------------------------

- When you are finished marking the cytosol, make the following changes

  - Select and hide the **Cytosol** [**h**].
  - Add a new boundary named **Mitochondria** and set the color to magenta.
  - Select one face on each mitochondria [**Shift+RMB**] and Select Linked
    [**Ctrl+l**]
  - Use **Assign** to assign the selected faces to be in the mitochondria.
  - When finished, hide the mitochondria [**h**] and proceed with marking the
    t-tubule (**TT**. Set color to blue) and sarcoplasmic reticulum (**SR**.
    Set color to yellow). We chose the two letter abbreviations because
    boundary names cannot contain special characters or spaces (underscores are
    OK).

  .. image:: ../images/all_marked.png

- CHECKPOINT: Save your progress to: **tt-sr-mit.all_marked.blend**

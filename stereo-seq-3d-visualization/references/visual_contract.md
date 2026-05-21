# Visual Contract

Use this checklist for 3D ST figure requests.

## Required Inputs

- `x`, `y`, and either `z` or `slice` with explicit `--slice-spacing`.
- A label, continuous value, or color column.
- Coordinate orientation: whether image-style y should be inverted.
- Slice order if slice names are not naturally sortable.

## Standard Outputs

- `*_view_elev*_azim*.png`: one static view per camera angle.
- `*_contact_sheet.png`: multi-angle contact sheet when more than one view is requested.
- `*_rotation.gif`: optional animation when requested.
- `*_single_slice_*.png`: optional black-background 2D slice panels.
- `*_provenance.json`: parameters, tool versions, input columns, and source-code provenance.

## View Presets

- `front`: elevation 12, azimuth -90.
- `side`: elevation 10, azimuth 0.
- `top`: elevation 90, azimuth -90.
- `oblique`: elevation 25, azimuth -60.
- `poster`: elevation 18, azimuth -45.

## Quality Checks

- Image is non-empty and points are visible on a black background.
- Point size does not turn dense regions into solid blocks unless that is intended.
- Slice spacing does not visually collapse sections or exaggerate anatomy without being reported.
- Categorical colors are stable across all snapshots and slices.
- Continuous colormap and normalization are identical across all views.
- Axes are hidden for presentation snapshots, but the report preserves coordinate orientation and camera angles.

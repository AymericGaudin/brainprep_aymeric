##########################################################################
# NSAp - Copyright (C) CEA, 2021 - 2026
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

"""
Synb0-DISCO functions.
"""

from ..decorators import (
    CoerceparamsHook,
    CommandLineWrapperHook,
    LogRuntimeHook,
    OutputdirHook,
    SignatureHook,
    step,
)
from ..typing import (
    Directory,
    File,
)
from ..utils import (
    bvecbval_from_file,
    sidecar_from_file,
)


@step(
    hooks=[
        CoerceparamsHook(),
        OutputdirHook(),
        LogRuntimeHook(
            bunched=False
        ),
        CommandLineWrapperHook(),
        SignatureHook(),
    ]
)
def synthb0(
        t1_file: File,
        dwi_file: File,
        workspace_dir: Directory,
        output_dir: Directory,
        entities: dict,
        mask_file: File | None = None
    ) -> tuple[list[list[str]], tuple[File]]:
    """
    Apply Synb0-DISCO.

    The Synb0-DISCO tool aims to enable susceptibility distortion correction
    with historical and/or limited datasets that do not include specific
    sequences for distortion correction (i.e. reverse phase-encoded scans). In
    short, the tool synthesizes an "undistorted" b=0 image that matches the
    geometry of structural T1w images and also matches the contrast from
    diffusion images. This synthesized 'undistorted' image can then be used
    in standard pipelines (i.e. TOPUP).

    Parameters
    ----------
    t1_file : File
        Path to the input T1w image file.
    dwi_file : File
        Path to the input diffusion weighted image file of one subject.
    workspace_dir: Directory
        Working directory with the workspace of the current processing.
    output_dir : Directory
        Directory where the generated images will be saved.
    entities : dict
        A dictionary of parsed BIDS entities including modality.
    mask_file: File | None
        Binary brain mask image file. If this parameter is not set, the mask
        is computed using bet.
        Default None.

    Returns
    -------
    commands : list[list[str]]
        Synb0-DISCO preprocessing command-lines.
    outputs : tuple[File]
        - b0_file : File - synthesized b=0 image as input to topup.
        - sidecar_file : File - associated JSON sidecar with reverse phase
          encoding and an effective echo spacing of 0 (infinite bandwidth).
    """
    basename = (
        f"sub-{entities['sub']}_"
        f"ses-{entities['ses']}_"
        f"run{entities['run']}"
    )

    fmap_dir = output_dir / "fmap"
    work_dir = workspace_dir / "work"
    for dir_ in (fmap_dir, work_dir, workspace_dir):
        dir_.mkdir(parents=True, exist_ok=True)

    sidecar_dwi_file = sidecar_from_file(dwi_file)
    _, bval_file = bvecbval_from_file(dwi_file)

    b0_file = fmap_dir / f"{basename}_dir-REVERSE_epi.nii.gz"
    sidecar_file = fmap_dir / f"{basename}_dir-REVERSE_epi.json"

    commands = [
        [
            "cp",
            str(t1_file),
            str(workspace_dir / "T1.nii.gz"),
        ] if mask_file is None else
        [
            "fslmaths",
            str(t1_file),
            "-mas", str(mask_file),
            str(workspace_dir / "T1.nii.gz"),
        ],
        [
            "fslroi",
            str(dwi_file),
            str(workspace_dir / "b0.nii"),
            (
                "$(awk "
                f"'{{for(i=1;i<=NF;i++) if($i<50) {{print i-1; exit}}}}' "
                f"{bval_file}"
                ")"
            ),
            "1",
        ],
        [
            "synb0",
            str(workspace_dir / "b0.nii"),
            str(workspace_dir / "T1.nii.gz"),
            str(work_dir),
            "mni_icbm152_t1_tal_nlin_asym_09c.nii.gz",
        ] if mask_file is None else
        [
            "synb0",
            str(workspace_dir / "b0.nii"),
            str(workspace_dir / "T1.nii"),
            str(work_dir),
            "mni_icbm152_t1_tal_nlin_asym_09c_mask.nii.gz",
        ],
        [
            "cp",
            str(work_dir / "b0_u.nii.gz"),
            str(b0_file),
        ],
        [
            "cp",
            str(sidecar_dwi_file),
            str(sidecar_file),
        ],
        [
            "sed",
            "-i",
            "-E",
            's/("TotalReadoutTime":[[:space:]]*)[0-9.]+/\\10/g',
            str(sidecar_file),
        ],
        [
            "sed",
            "-i",
            "-E",
            (
                's/("PhaseEncodingDirection":[[:space:]]*"[^"]+)-"/\\1"/g; '
                't; s/("PhaseEncodingDirection":[[:space:]]*"[^"-]+")/\\1-"/g'
            ),
            str(sidecar_file),
        ],
    ]

    return commands, [
        b0_file,
        sidecar_file,
    ]

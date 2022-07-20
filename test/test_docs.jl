# This file is a part of AffineTransformations.jl, licensed under the MIT License (MIT).

using Test
using AffineTransformations
import Documenter

Documenter.DocMeta.setdocmeta!(
    AffineTransformations,
    :DocTestSetup,
    :(using AffineTransformations);
    recursive=true,
)
Documenter.doctest(AffineTransformations)

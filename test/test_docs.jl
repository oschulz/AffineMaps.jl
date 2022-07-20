# This file is a part of AffineMaps.jl, licensed under the MIT License (MIT).

using Test
using AffineMaps
import Documenter

Documenter.DocMeta.setdocmeta!(
    AffineMaps,
    :DocTestSetup,
    :(using AffineMaps);
    recursive=true,
)
Documenter.doctest(AffineMaps)

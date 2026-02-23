{
  description = "hmm - heightmap meshing utility";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
      in
      {
        packages = {
          default = pkgs.stdenv.mkDerivation {
            pname = "hmm";
            version = "0.1.0";

            src = ./.;

            nativeBuildInputs = with pkgs; [
              gnumake
            ];

            buildInputs = with pkgs; [
              glm
            ];

            makeFlags = [ "C=${pkgs.stdenv.cc}/bin/c++" ];

            buildPhase = ''
              runHook preBuild
              make release
              runHook postBuild
            '';

            installPhase = ''
              runHook preInstall
              mkdir -p $out/bin
              cp bin/release/hmm $out/bin/
              runHook postInstall
            '';

            meta = with pkgs.lib; {
              description = "Heightmap meshing utility";
              longDescription = ''
                hmm is a heightmap meshing utility that converts grayscale heightmap
                images into 3D meshes using a fast polygonal approximation algorithm.
                The meshes satisfy the Delaunay condition and can satisfy a specified
                maximal error or maximal number of triangles or vertices.
              '';
              homepage = "https://github.com/fogleman/hmm";
              license = licenses.mit;
              platforms = platforms.unix;
              mainProgram = "hmm";
            };
          };
        };

        devShells.default = pkgs.mkShell {
          buildInputs = with pkgs; [
            glm
            gnumake
            gcc
          ];

          shellHook = ''
            echo "hmm development environment"
            echo "Run 'make' to build the project"
          '';
        };
      }
    );
}

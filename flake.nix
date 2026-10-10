{
  outputs = inputs:
    inputs.flake-parts.lib.mkFlake {inherit inputs;} ({...}: {
      systems = ["x86_64-linux" "aarch64-linux"];

      perSystem = {pkgs, ...}: {
        devShells.default = pkgs.mkShell {
          packages = with pkgs; [
            python314
            ruff
            uv
          ];

          env = {
            UV_PYTHON_PREFERENCE = "only-system";
            UV_PYTHON_DOWNLOADS = "never";
            UV_PROJECT_ENVIRONMENT = ".venv";

            # stdenv adds `-rpath $PWD/outputs/out/lib` to NIX_LDFLAGS, which the linker wrapper 
            # splits at spaces: building rainbow-api's C extension then fails in a checkout
            # whose path contains spaces.
            NIX_NO_SELF_RPATH = "1";

            # PyPI wheels (numpy, scipy) link against these.
            LD_LIBRARY_PATH = pkgs.lib.makeLibraryPath (with pkgs; [
              stdenv.cc.cc.lib # libstdc++, libgcc_s
              zlib
            ]);
          };
        };
      };
    });

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-parts.url = "github:hercules-ci/flake-parts";
  };
}

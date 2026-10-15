{
  description = "PacIFiC (source-tree)";

  # Optional: helps CI/users if they haven't enabled flakes globally
  nixConfig = {
    experimental-features = [
      "nix-command"
      "flakes"
    ];
  };

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

  outputs =
    { self, nixpkgs }:
    let
      systems = [
        "x86_64-linux"
        "aarch64-linux"
      ];
      forAllSystems = f: nixpkgs.lib.genAttrs systems (system: f system);
    in
    {
      overlays.default =
        final: prev:
        let
          cudaToolkit =
            if final.stdenv.isLinux && final ? cudaPackages && final.cudaPackages ? cudatoolkit then
              final.cudaPackages.cudatoolkit
            else
              null;
          cudaNvcc =
            if final.stdenv.isLinux && final ? cudaPackages && final.cudaPackages ? cuda_nvcc then
              final.cudaPackages.cuda_nvcc
            else
              null;
        in
        {
          pacific = final.callPackage ./pacific.nix {
            inherit cudaToolkit cudaNvcc;
          };
        };

      packages = forAllSystems (
        system:
        let
          pkgs = import nixpkgs {
            inherit system;
            overlays = [ self.overlays.default ];
            config.allowUnfree = true;
          };
        in
        {
          pacific = pkgs.pacific;
          default = pkgs.pacific;
        }
      );

      devShells = forAllSystems (
        system:
        let
          pkgs = import nixpkgs {
            inherit system;
            overlays = [ self.overlays.default ];
            config.allowUnfree = true;
          };
          cudaToolkit =
            if pkgs.stdenv.isLinux && pkgs ? cudaPackages && pkgs.cudaPackages ? cudatoolkit then
              pkgs.cudaPackages.cudatoolkit
            else
              null;
          cudaNvcc =
            if pkgs.stdenv.isLinux && pkgs ? cudaPackages && pkgs.cudaPackages ? cuda_nvcc then
              pkgs.cudaPackages.cuda_nvcc
            else
              null;
        in
        {
          default = pkgs.mkShell {
            inputsFrom = [ pkgs.pacific ];
            packages =
              (with pkgs; [
                gdb
                petsc
                doxygen
                python313Packages.mkdocs
                apptainer
              ])
              ++ pkgs.lib.optionals (cudaToolkit != null) [ cudaToolkit ]
              ++ pkgs.lib.optionals (cudaNvcc != null) [ cudaNvcc ];
            shellHook = ''
              export PETSC_DIR="${pkgs.petsc}"
              export HYPRE_DIR="${pkgs.hypre}"
            ''
            + pkgs.lib.optionalString (cudaToolkit != null) ''
              export CUDA_PATH="${cudaToolkit}"
              export CUDA_HOME="${cudaToolkit}"
            ''
            + pkgs.lib.optionalString (cudaNvcc != null) ''
              export CUDACXX="${cudaNvcc}/bin/nvcc"
            '';
          };
        }
      );

      # Nice: `nix flake check` builds it
      checks = forAllSystems (system: {
        pacific = self.packages.${system}.pacific;
      });
    };
}

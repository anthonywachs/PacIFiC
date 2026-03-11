{
  description = "PacIFiC (source-tree)";

  # Optional: helps CI/users if they haven't enabled flakes globally
  nixConfig = {
    experimental-features = [ "nix-command" "flakes" ];
  };

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

  outputs = { self, nixpkgs }:
    let
      systems = [ "x86_64-linux" "aarch64-linux" ];
      forAllSystems = f: nixpkgs.lib.genAttrs systems (system: f system);
    in
    {
      overlays.default = final: prev: {
        pacific = final.callPackage ./pacific.nix { };
      };

      packages = forAllSystems (system:
        let
          pkgs = import nixpkgs {
            inherit system;
            overlays = [ self.overlays.default ];
          };
        in
        {
          pacific = pkgs.pacific;
          default = pkgs.pacific;
        });

      devShells = forAllSystems (system:
        let
          pkgs = import nixpkgs {
            inherit system;
            overlays = [ self.overlays.default ];
          };
        in
        {
          default = pkgs.mkShell {
            inputsFrom = [ pkgs.pacific ];
            packages = with pkgs; [
              petsc
              hypre
              gdb
              doxygen
              python313Packages.mkdocs
              apptainer
            ];
            shellHook = ''
              export PETSC_DIR="${pkgs.petsc}"
              export HYPRE_DIR="${pkgs.hypre}"
            '';
          };
        });

      # Nice: `nix flake check` builds it
      checks = forAllSystems (system: {
        pacific = self.packages.${system}.pacific;
      });
    };
}
